import math
import numpy as np

import argparse
import os
from itertools import product

def round_prob(p: float, num_bit: int):
    bits = 0
    list_bits = []
    for i in range(num_bit):
        p *= 2
        if p >= 1:
            p -= 1
            bits |= 1 << (num_bit - i - 1)
            list_bits.append(1)
        else:
            list_bits.append(0)        
        # remove the 0s at the end
        while len(list_bits) > 0 and list_bits[-1] == 0:
            list_bits.pop()
    return bits / (1 << num_bit), len(list_bits)

def bc23_get_cost(eps: float, delta: float):
    # compute the number of bits
    nb = math.ceil(100 * math.log(2/delta) / (eps ** 2))
    l1 = math.sqrt(nb/2) / math.sqrt(math.pi)
    return eps, delta, nb, l1

def safe_geom_to_ber(p_geom: float, i: int):
    try:
        return 1 / (1 + math.pow(1 - p_geom, -1 << i))
    except OverflowError:
        return 0

def vdlm_get_priv_cost(t: float, log_range: int, prec: int):
    # probability of returning 0
    pz_ = (math.exp(1/t) - 1) / (math.exp(1/t) + 1)
    pz, num_bit = round_prob(pz_, prec)
    qz = 1 - pz
    num_bit += 1 # the one bit for the sign

    p_geom = 1 - math.exp(-1/t)
    p_list_ = [safe_geom_to_ber(p_geom, i) for i in range(log_range)]
    p_list = []
    for p_ in p_list_:
        p, temp_num_bits = round_prob(p_, prec)
        p_list.append(p)
        num_bit += temp_num_bits
    # remove the terms with p = 0
    while p_list[-1] == 0:
        p_list.pop()
        log_range -= 1
    q_list = [1 - p for p in p_list]
    log_p_list = [math.log(p) for p in p_list]
    log_q_list = [math.log(1-p) for p in p_list]

    # cumsum of log_p and log_q
    cumsum_log_p = [0] + np.cumsum(log_p_list).tolist()
    cumsum_log_q = [0] + np.cumsum(log_q_list).tolist()

    # compute epsilon and delta 
    eps = math.log(pz) + math.log(2) - math.log(qz) - cumsum_log_q[-1]
    for i in range(log_range):
        eps_temp = log_q_list[i] - log_p_list[i] + cumsum_log_p[i] - cumsum_log_q[i]
        eps = max(eps, eps_temp)
    delta = (1 - pz) * math.exp(cumsum_log_p[-1]) / 2 

    l1_expectation = 0
    for i in range(log_range):
        l1_expectation += (1 << i) * p_list[i]
    l1_expectation = qz * (l1_expectation + 1)
    

    return eps, delta, num_bit, l1_expectation

argparser = argparse.ArgumentParser()
argparser.add_argument('--eps', type=float, required=True)
argparser.add_argument('--delta', type=float, required=True)
argparser.add_argument('--count_search_space', type=int, required=True)
argparser.add_argument('--prec_search_space', type=int, required=True)
argparser.add_argument('--eps_tolerance', type=float, default=1.1)

if __name__ == "__main__":
    args = argparser.parse_args()
    _, _, bc23_num_bit, bc23_l1 = bc23_get_cost(args.eps, args.delta)

    if not os.path.isdir('logs'):
        os.mkdir('logs')
        
    if not os.path.isfile('logs/vddlm-anal.csv'):
        with open('logs/vddlm-anal.csv', 'w') as f:
            f.write('algorithm,eps,delta,num_coin,l1\n')

    with open('logs/vddlm-anal.csv', 'a') as f:
        f.write(f'bc23,{args.eps},{args.delta},{bc23_num_bit},{bc23_l1}\n')

    
    for log_range, prec in product(range(1, args.count_search_space), range(1, args.prec_search_space)):
        try:
            eps, delta, num_bit, l1_expectation = vdlm_get_priv_cost(1/args.eps, log_range, prec)
            if eps < args.eps * args.eps_tolerance and delta < args.delta:
                with open('logs/vddlm-anal.csv', 'a') as f:
                    f.write(f'ours,{eps},{delta},{num_bit},{l1_expectation}\n')
                break
        except:
            pass
