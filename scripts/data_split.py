import shutil
import os
import sys
import numpy as np
from concurrent.futures import ThreadPoolExecutor

path = sys.argv[1]

system = {
    "NaCl": {"Na", "Cl"},
    "NaK": {"Na", "K"},
    "NaClK": {"Na", "Cl", "K"},
    "KCl": {"K", "Cl"},
    "Na": {"Na"},
    "Cl": {"Cl"},
    "K": {"K"},
    "others": {"none"}
}

def sites(lines):
    number_of_atoms_cell = [line for line in lines if 'number of atoms/cell' in line] 
    Nat = int(number_of_atoms_cell[0].split()[4]) if number_of_atoms_cell else None
    keyword = ''
    Type = [] 
    Line_count = 0       
    for line in lines:
        if 'site n.' in line:
            keyword = 'site n.'
            continue
        if keyword == 'site n.':
            tokens = line.split()
            if len(tokens) > 1:
                Type.append(tokens[1])
            Line_count += 1
            if Line_count == Nat:
                keyword = ''
                Line_count = 0
    return Type

def move_to_dir(path):
    for dirname in system.keys():
        dir_path = os.path.join(path, dirname)
        if not os.path.exists(dir_path):
            print(f"Creating directory: {dir_path}")
            os.makedirs(dir_path)
        else:
            print(f"Directory already exists: {dir_path}")

def move_file(path, dirname, filename):        
    dir_path = os.path.join(path, dirname)      
    src_path = os.path.join(path, filename)
    dest_path = os.path.join(dir_path, filename)
    shutil.copy(src_path, dest_path)

def process_file(path, filename):
    if filename.endswith('.out') or 'out' in filename:
        full_path = os.path.join(path, filename)
        with open(full_path, 'r') as f:
            lines = f.readlines()
        unique_elements = np.unique(sites(lines))
        elements_set = set(unique_elements)
        matching_keys = [key for key, value in system.items() if value == elements_set]
        if matching_keys:
            dirname = matching_keys[0]
            move_file(path, dirname, filename)
            
        else:
            dirname = "others" 
            move_file(path, dirname, filename)


def read_files(path):
    move_to_dir(path)
    filenames = os.listdir(path)
    with ThreadPoolExecutor() as executor:
        for filename in filenames:
            executor.submit(process_file, path, filename)

read_files(path)
