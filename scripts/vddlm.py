import os, sys
import math
import glog

import argparse
import subprocess
import numpy as np

import tempfile
from utils import save_int, load_int
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
    while len(list_bits) > 0 and list_bits[-1] == 0:
        list_bits.pop()
    return bits / (1 << num_bit), list_bits

def bc23_get_cost(eps: float, delta: float):
    # compute the number of bits
    nb = math.ceil(10 * math.log(2/delta) / (eps ** 2))
    l1 = math.sqrt(nb/2) / math.sqrt(math.pi)
    l2 = math.sqrt(nb) / 2
    return eps, delta, nb, l1, l2

def safe_geom_to_ber(p_geom: float, i: int):
    try:
        return 1 / (1 + math.pow(1 - p_geom, -1 << i))
    except OverflowError:
        return 0
    
def vdlm_get_priv_cost(t: float, log_range: int, prec: int):
    prob_bits = {}
    # probability of returning 0
    pz_ = (math.exp(1/t) - 1) / (math.exp(1/t) + 1)
    pz, prob_bits['z'] = round_prob(pz_, prec)
    qz = 1 - pz

    p_geom = 1 - math.exp(-1/t)
    p_list_ = [safe_geom_to_ber(p_geom, i) for i in range(log_range)]
    p_list = [None] * len(p_list_)
    for i, p_ in enumerate(p_list_):
        p_list[i], prob_bits[i] = round_prob(p_, prec)
    # remove the terms with p = 0
    while p_list[-1] == 0:
        p_list.pop()
        log_range -= 1
        prob_bits.pop(len(p_list))
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

    g_bitwise_exp = np.array([(1 << i) * p_list[i] for i in range(log_range)])
    g_sq_bitwise_exp = g_bitwise_exp[:, None] * g_bitwise_exp[None, :]

    # exp of geom
    g_exp = g_bitwise_exp.sum()
    g_sq_exp = g_sq_bitwise_exp.sum()

    # compute the expectation of L1
    l1 = qz * (g_exp + 1)
    l2 = math.sqrt(qz * (g_sq_exp + 2 * g_exp + 1))

    num_bit = 1 + sum([len(_) for _ in prob_bits.values()])

    return eps, delta, num_bit, l1, l2, prob_bits

parser = argparse.ArgumentParser()
parser.add_argument('--dim', type=int) # not required in anal mode
parser.add_argument('--eps', type=float, required=True)
parser.add_argument('--delta', type=float, required=True)
parser.add_argument('--count_range', type=int) # not required in anal mode
parser.add_argument('--log_noise_range_search_space', type=int, required=True)
parser.add_argument('--prec_search_space', type=int, required=True)
parser.add_argument('--eps_tolerance', type=float, default=1.05)
parser.add_argument('--anal', action = 'store_true')
parser.add_argument('--nser', type=int, default = 2)

if __name__ == '__main__':
    args = parser.parse_args()

    if args.anal:
        _, _, bc23_num_bit, bc23_l1, bc23_l2 = bc23_get_cost(args.eps, args.delta)
        if not os.path.isdir('logs'):
            os.mkdir('logs')
            
        if not os.path.isfile('logs/vddlm-anal.csv'):
            with open('logs/vddlm-anal.csv', 'w') as f:
                f.write('algorithm,eps,delta,num_coin,l1,l2\n')

        with open('logs/vddlm-anal.csv', 'a') as f:
            f.write(f'bc23,{args.eps},{args.delta},{bc23_num_bit},{bc23_l1},{bc23_l2}\n')

    eps, delta, num_bit, l1, l2, prob_bits = None, None, None, None, None, None
    for log_noise_range, prec in product(range(1, args.log_noise_range_search_space), range(1, args.prec_search_space)):
        try:
            eps, delta, num_bit, l1, l2, prob_bits = vdlm_get_priv_cost(1/args.eps, log_noise_range, prec)
            if eps < args.eps * args.eps_tolerance and delta < args.delta:
                if args.anal:
                    with open('logs/vddlm-anal.csv', 'a') as f:
                        f.write(f'ours,{eps},{delta},{num_bit},{l1},{l2}\n')
                break
        except Exception as e:
            if isinstance(e, IndexError) or isinstance(e, ValueError):
                pass
            else:
                raise e

    if not args.anal:
    
        # get a temporary directory
        with tempfile.TemporaryDirectory() as tempdir:
            with open(os.path.join(tempdir, 'prob_bits.txt'), 'w') as f:
                f.write('z ' + ''.join([str(_) for _ in prob_bits['z']]) + '\n')
                for i in range(len(prob_bits) - 1):
                    f.write(f'{i} ' + ''.join([str(_) for _ in prob_bits[i]]) + '\n')

            # Generate count data in range [0, count_range]
            count_data = np.random.randint(0, args.count_range, args.dim)
            count_data = count_data.astype(np.int32)
            count_data_path = os.path.join(tempdir, 'count_data.bin')
            count_noisy_data_path = os.path.join(tempdir, 'count_noisy_data.bin')
            count_data.tofile(count_data_path)

            output = subprocess.check_output(f'./build/src/vddlm {args.nser} {args.dim} {count_data_path} {count_noisy_data_path} {tempdir}/prob_bits.txt', shell=True)
            output_str = output.decode('utf-8')

            # get the noisy count data
            noise_count_data = load_int(count_noisy_data_path)

            # get average L1 and L2 error between the original and noisy count data
            l1_diff = np.abs(count_data - noise_count_data).sum() 
            l2_diff = np.sqrt(np.square(count_data - noise_count_data).sum()) 

            if not os.path.isdir('logs'):
                os.mkdir('logs')
            
            if not os.path.isfile('logs/vddlm.csv'):
                with open('logs/vddlm.csv', 'w') as f:
                    f.write('nser,dim,eps,delta,num_bit,setup,computing,proving,verifying,l1,l2\n')
            
            log_str = f'{args.nser},{args.dim},{eps},{delta},{num_bit},'
            rts = [line.split(' ')[-2] for line in output_str.split('\n') if line]
            assert(len(rts) == 4)
            log_str += ','.join(rts)
            log_str += f',{l1_diff},{l2_diff}'
            with open('logs/vddlm.csv', 'a') as f:
                f.write(log_str + '\n')



