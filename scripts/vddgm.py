import os, sys
import math
import glog

import argparse
import subprocess
import numpy as np

import tempfile
from utils import save_int, load_int
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument('--nser', type=int, default = 2)
parser.add_argument('--dim', type=int, required = True) 
parser.add_argument('--count_range', type=int, default=114514)
parser.add_argument('--sigma', type=float)
parser.add_argument('--prob_nbit', type=int, default=20)

if __name__ == '__main__':
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tempdir:

        count_data = np.random.randint(0, args.count_range, args.dim)
        count_data = count_data.astype(np.int32)
        count_data_path = os.path.join(tempdir, 'count_data.bin')
        count_noisy_data_path = os.path.join(tempdir, 'count_noisy_data.bin')
        count_data.tofile(count_data_path)

        command = f'./build/src/vddgm {args.nser} {args.dim} {args.sigma} {args.prob_nbit} {count_data_path} {count_noisy_data_path}'
        output_str = subprocess.check_output(command, shell=True).decode('utf-8')

        # get the noisy count data
        noise_count_data = load_int(count_noisy_data_path)

        # get average L1 and L2 error between the original and noisy count data
        l1_diff = np.abs(count_data - noise_count_data).sum() 
        l2_diff = np.sqrt(np.square(count_data - noise_count_data).sum()) 

        if not os.path.isdir('logs'):
            os.mkdir('logs')
        
        if not os.path.isfile('logs/vddgm.csv'):
            with open('logs/vddgm.csv', 'w') as f:
                f.write('nser,dim,sigma,prob_nbit,setup,computing,proving,verifying,comm,l1,l2\n')

        log_str = f'{args.nser},{args.dim},{args.sigma},{args.prob_nbit},'
        rts = [line.split(' ')[-2] for line in output_str.split('\n') if line]
        assert(len(rts) == 5)
        log_str += ','.join(rts)
        log_str += f',{l1_diff},{l2_diff}'
        with open('logs/vddgm.csv', 'a') as f:
            f.write(log_str + '\n')