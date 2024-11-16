import os, sys
import math
import glog

import argparse
import subprocess
import numpy as np

import tempfile
from utils import save_int, load_int

parser = argparse.ArgumentParser()
parser.add_argument('--dim', type=int, required=True)
parser.add_argument('--eps', type=float, required=True)
parser.add_argument('--count_range', type=int, required=True)
parser.add_argument('--noise_log_range', type=int, required=True)
parser.add_argument('--prec', type=int, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    
    # get a temporary directory
    with tempfile.TemporaryDirectory() as tempdir:

        # Generate count data in range [0, count_range]
        count_data = np.random.randint(0, args.count_range, args.dim)
        count_data = count_data.astype(np.int32)
        count_data_path = os.path.join(tempdir, 'count_data.bin')
        count_noisy_data_path = os.path.join(tempdir, 'count_noisy_data.bin')
        count_data.tofile(count_data_path)

        output = subprocess.check_output(f'./build/src/vddlm {count_data_path} {count_noisy_data_path} {args.dim} {args.eps} {args.noise_log_range} {args.prec}', shell=True)
        output_str = output.decode('utf-8')

        # get the noisy count data
        noise_count_data = load_int(count_noisy_data_path)

        # get average L1 and L2 error between the original and noisy count data
        l1_diff = np.abs(count_data - noise_count_data).sum() / args.dim
        l2_diff = np.sqrt(np.square(count_data - noise_count_data).sum()) / args.dim

        if not os.path.isdir('logs'):
            os.mkdir('logs')
        
        if not os.path.isfile('logs/vddlm.csv'):
            with open('logs/vddlm.csv', 'w') as f:
                f.write('dim,eps,noise_log_range,prec,setup,computing,proving,verifying,l1,l2\n')
        
        log_str = f'{args.dim},{args.eps},{args.noise_log_range},{args.prec},'
        rts = [line.split(' ')[-2] for line in output_str.split('\n') if line]
        assert(len(rts) == 4)
        log_str += ','.join(rts)
        log_str += f',{l1_diff},{l2_diff}'
        with open('logs/vddlm.csv', 'a') as f:
            f.write(log_str + '\n')



