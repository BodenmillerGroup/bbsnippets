#!/usr/bin/env python3
import h5py

h5_file = '/path/to/file.h5'

def print_path(prefix, name, obj):
    if isinstance(obj, h5py.Dataset):
        print(prefix + name, obj.shape)
    else:
        for sub_name in obj:
            print_path(prefix + name + '/', sub_name, obj[sub_name])

if __name__ == '__main__':
    with h5py.File(h5_file, 'r') as f:
        for name in f:
            print_path('/', name, f[name])
