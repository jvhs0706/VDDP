import os, sys
import math
import glog

import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--eps', type=float, required=True)
parser.add_argument('--log_prec', type=int, required=True)
parser.add_argument('--num_class', type=int, required=True)

def get_vrr_params(eps: float, log_prec: int, num_class: int):
    assert log_prec > 0 and num_class > 1
    Omega = 1 << log_prec
    B_ = Omega / (math.exp(eps) + num_class - 1)
    B = round(B_)
    assert isinstance(B, int)
    A = Omega - (num_class - 1) * B
    return Omega, A, B, math.log(A / B)
    

if __name__ == '__main__':
    args = parser.parse_args()
    try:
        Omega, A, B, actual_eps = get_vrr_params(args.eps, args.log_prec, args.num_class)
    except ZeroDivisionError:
        glog.error(f'Invalid parameters with eps = {args.eps}, log_prec = {args.log_prec}, num_class = {args.num_class}!')
        exit(1)

    output = subprocess.check_output(f'./build/src/vrr {args.num_class} {Omega} {A} {B}', shell=True)
    output_str = output.decode('utf-8')

    if not os.path.isdir('logs'):
        os.mkdir('logs')
    
    if not os.path.isfile('logs/vrr.csv'):
        with open('logs/vrr.csv', 'w') as f:
            f.write('eps,log_prec,num_class,actual_eps,setup,committing,computing,proving,verifying,comm\n')

    
    log_str = f'{args.eps},{args.log_prec},{args.num_class},{actual_eps},'
    rts = [line.split(' ')[-2] for line in output_str.split('\n') if line]
    assert(len(rts) == 6)
    log_str += ','.join(rts)
    with open('logs/vrr.csv', 'a') as f:
        f.write(log_str + '\n')


