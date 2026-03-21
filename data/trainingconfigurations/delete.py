def delete():
    import sys
    import os
    import csv
    import numpy as np
    import numpy
    import shutil
    import glob
    import math
    import subprocess

    qe_input_dir ='/home/karimzongo/Desktop/T6'
    Ry = 13.605693012183622
    Ang = 0.5291772105638411
    F_conversion = Ry/Ang
    S_conversion = Ry/(Ang**3)
    Total = 0

#    files = glob.glob(r'./*')
#    for items in files:
       # os.remove(items)
    # Unity for conversion to eV or angstrom
    for filename in os.listdir(qe_input_dir):
        if filename.startswith('out.'):
            Total = Total + 1
            with open(os.path.join(qe_input_dir, filename), 'r') as txtfile:
                LINES = txtfile.readlines()
                #print(lines)
                E = []
                Positions = []
                Prim_vector = []
                CELL = []
                lattice = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float)
                lines = []

                index = 0
                for line in LINES:
                    line = line.strip()
                    if line:
                            lines.append(line)
                    if line == 'JOB DONE.':
                        lines.append(line)
                        break
                l_file = len(lines)
                l_job_done = l_file - 2

                #print(lines[l_file - 2])
                for line in lines:
                    if ('convergence NOT achieved after ' in line):
                        sys.stdout.write('rm ')
                        print(' ', filename)
                #elif (lines[l_job_done] == 'JOB DONE.'):
                if (lines[l_job_done] != 'JOB DONE.'):
                    sys.stdout.write('rm ')
                    print(' ', filename)
        txtfile.close()

                        #print(Positions_P)
                        #print(Cell_dim)
                        #print(lattice)
                        #print(Energy)
                        #print(Forces_F)
                        #print(Stress_S)
                        #print(Energy)

delete()

