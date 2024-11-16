import numpy as np

def save_int(arr: np.ndarray, path):
    if path[-4:] != '.bin':
        raise ValueError('Path must end with .bin')
    assert arr.dtype == np.int32, 'Array must be of type np.int32'
    arr.tofile(path)

def save_long(arr: np.ndarray, path):
    if path[-4:] != '.bin':
        raise ValueError('Path must end with .bin')
    assert arr.dtype == np.int64, 'Array must be of type np.int64'
    arr.tofile(path)

def load_int(path):
    if path[-4:] != '.bin':
        raise ValueError('Path must end with .bin')
    return np.fromfile(path, dtype=np.int32)

def load_long(path):
    if path[-4:] != '.bin':
        raise ValueError('Path must end with .bin')
    return np.fromfile(path, dtype=np.int64)