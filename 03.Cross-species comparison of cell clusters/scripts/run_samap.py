#!/usr/bin/env python3

#import os
#from samalg import SAM

from scanpy import read_h5ad
from pathlib import Path
import pickle
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--abbr1', type=str, required=True, help='Abbreviation of 1st species')
    parser.add_argument('-A', '--abbr2', type=str, required=True, help='Abbreviation of 2nd species')
    parser.add_argument('-k', '--key1', type=str, required=True, help='Key column name in metadata of 1st scanpy object')
    parser.add_argument('-K', '--key2', type=str, required=True, help='Key column name in metadata of 2nd scanpy object')
    parser.add_argument('-f', '--file1', type=str, required=True, help='H5ad file of 1st scanpy object')
    parser.add_argument('-F', '--file2', type=str, required=True, help='H5ad file of 2st scanpy object')
    parser.add_argument('-b', '--blast-out', type=str, required=True, help='Directory stored blast output files')
    parser.add_argument('-o', '--output-dir', type=str, required=True, help='Output directory')
    return parser.parse_args()


def check_file_exists(filenames):
    for k,f in filenames.items():
        if Path(f).exists() is False:
            raise FileNotFoundError(f'Error: H5ad file: {f} not found for {k}.\n')


def run_samap(id1, id2, keys, filenames, blast_out_dir, out_dir):
    def check_keys(file, key, out_dir):
        new_file = str(Path(out_dir) / Path(file).name)
        obj = read_h5ad(file)
        if key not in obj.obs.columns:
            raise Exception(f'Key {key} cannot be found in {file}.')
        obj.obs[key] = obj.obs[key].astype('str').astype('category')
        obj.write_h5ad(new_file)
        return new_file
        
    for i in [id1, id2]:
        filenames[i] = check_keys(filenames[i], keys[i], out_dir)
    out_pkl = Path(f'{out_dir}/samap.pkl')
    if out_pkl.exists():
        print(f'SAMap result object already existed. Running SAMap have been skipped.\n    {out_pkl}')
        return
    print(f'SAMap result object did not exist. Continue to run SAMap...')
    from samap.mapping import SAMAP
    sm = SAMAP(
            filenames,
            f_maps = blast_out_dir,
            keys = keys,
            save_processed = False #if False, do not save the processed results to `*_pr.h5ad`
        )
    sm.run(neigh_from_keys = {id1:True, id2:True})
    with open(f'{out_dir}/samap.pkl', 'wb') as f:
        pickle.dump(sm, f)


def main():
    import argparse
    args = parse_arguments()
    id1 = args.abbr1
    id2 = args.abbr2
    out_dir = args.output_dir
    Path(out_dir).mkdir(exist_ok=True)
    keys = {id1: args.key1, id2: args.key2}
    print('Keys:', keys)
    filenames = {id1: args.file1, id2: args.file2}
    print('Files:', filenames)
    check_file_exists(filenames)

    blast_out_dir = args.blast_out
    blast_out_dir = blast_out_dir.rstrip('/') + '/'
    run_samap(id1, id2, keys, filenames, blast_out_dir, out_dir)


if __name__ == '__main__':
    main()
