import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
import random
import subprocess

filename = "mlip_input.cfg"
path=sys.argv[1]

def cluster(energy):

    Energies = np.array(energy)         
    X = Energies.reshape(-1,1)  
    db = DBSCAN(eps=0.2, min_samples=2).fit(X)

    labels = db.labels_
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print("Clusters:", n_clusters)
    
    return n_clusters, labels

def print_clusters(energies, labels):
    """
    Print each cluster label with its corresponding energy values.
    Noise (label -1) is listed separately.
    """
    energies = np.array(energies)
    labels = np.array(labels)
    unique_labels = sorted(set(labels))
    
    clean_energy = []
    close_energy = []
    
    for lbl in unique_labels:
        mask = (labels == lbl)
        cluster_vals = energies[mask]
        if lbl == -1:
            clean_energy.append(cluster_vals)
            print(f" Isolated points (configurations), count {len(cluster_vals)}:")
        else:
            close_energy.append(cluster_vals)
            print(f"Cluster {lbl}, count {len(cluster_vals)} (group of similar configurations):")
        print(cluster_vals)
        #print("---")

    flat_clean = [val for arr in clean_energy for val in arr]
    flat_close = [arr.tolist() for arr in close_energy]

    return flat_clean, flat_close
   

def extract_data_from_cfg(lines,disperse_energy, close_energy):

    current_block = []
    extracted_blocks = []

    for line in lines:
        stripped = line.strip()
        current_block.append(stripped)

        if stripped == "END_CFG":
            block = current_block
            current_block = []

            # Now extract data from the block:
            cfg = {"energy": None,"Supercell": [], "size": None, "stress": None, "atoms": []}

            for i, l in enumerate(block):
                if l == "Size":
                    cfg["size"] = int(block[i + 1])

                elif l == "Energy":
                    cfg["energy"] = float(block[i + 1])
                    
                elif l == "Supercell":
                
                     cell_infos = block[i + 1:]
                     
                     for vector in cell_infos:
                        
                         if not cell_infos or vector.startswith("AtomData"):
                         
                             break
                             
                         components = vector.split()
                         cell_vectors = {
                            "V1": float(components[0]),
                            "V2": float(components[1]),
                            "V3": float(components[2]),
                         }
                         cfg["Supercell"].append(cell_vectors)
                     
                elif l.startswith("PlusStress"):
                    stress_vals = block[i + 1].split()
                    cfg["stress"] = list(map(float, stress_vals))

                elif l.startswith("AtomData"):
                    # Parse atom data until a blank line or "Energy"
                    atom_lines = block[i + 1:]
                    for atom_line in atom_lines:
                        if not atom_line or atom_line.startswith("Energy"):
                            break
                        parts = atom_line.split()
                        atom = {
                            "id": int(parts[0]),
                            "type": int(parts[1]),
                            "x": float(parts[2]),
                            "y": float(parts[3]),
                            "z": float(parts[4]),
                            "fx": float(parts[5]),
                            "fy": float(parts[6]),
                            "fz": float(parts[7]),
                        }
                        cfg["atoms"].append(atom)

            extracted_blocks.append(cfg)    

    return extracted_blocks
        #rounded_disperse = [round(float(val), 6) for arr in disperse_energy for val in arr]
            
def cfg_energy(path):
   
    energies=[]
    keyword = ''
    with open(os.path.join(path, filename), 'r') as f:
    
       lines = f.readlines()
       
       for line in lines:
       
           if 'Energy' in line:
           
               keyword = 'Energy'
               continue
               
           if keyword == 'Energy':
           
               for x in line.split():
                
                   energies.append(float(x))
                   
                   keyword = ''
    
    return  lines, energies     
                   
def write_conf(size,supercell,atom_infos,Stress,Energy):
        
        lattice = []
        
        for vector in supercell:
           
           lattice.append(vector['V1'])
           lattice.append(vector['V2'])
           lattice.append(vector['V3'])        
        
        xx = Stress[0]
        yy = Stress[1]
        zz = Stress[2]
        yz = Stress[3]
        xz = Stress[4]
        xy = Stress[5]
               
        mlip_conf.write("BEGIN_CFG\n")
        mlip_conf.write(" Size\n")
        mlip_conf.write("    " + str(size) + "\n")
        mlip_conf.write(" Supercell\n")
        mlip_conf.write("    " + str(lattice[0]) + "     " + str(lattice[1]) + "      " + str(lattice[2]) + " \n")
        mlip_conf.write("    " + str(lattice[3]) + "     " + str(lattice[4]) + "      " + str(lattice[5]) + " \n")
        mlip_conf.write("    " + str(lattice[6]) + "     " + str(lattice[7]) + "      " + str(lattice[8]) + " \n")

        mlip_conf.write(" AtomData:  id type    cartes_x    cartes_y    cartes_z    fx    fy    fz\n")

        for i in range(0, size):
           
            atom = atom_infos[i]
            
            pos_x = atom['x']
            pos_y = atom['y']
            pos_z = atom['z']
            
            f_x = atom['fx']
            f_y = atom['fy']
            f_z = atom['fz']
            
            typ = atom['type']
            
            mlip_conf.write(
                "    " + str(i + 1) + "    " + str(typ) + "   " + str(pos_x) + "    " + str(pos_y) + "    " + str(
                    pos_z) + "    " + str(f_x) + "    " + str(f_y) + "    " + str(f_z) + "\n")
        
        mlip_conf.write(" Energy\n")
        mlip_conf.write("    " + str(Energy) + "\n")
        mlip_conf.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        mlip_conf.write(
            "    " + str(xx) + "   " + str(yy) + "    " + str(zz) + "    " + str(
                yz) + "    " + str(xz) + "    " + str(xy) + "\n")
        mlip_conf.write(" Feature EFS_by Qe\n")
        mlip_conf.write("END_CFG\n")
        mlip_conf.write("\n")

mlip_conf = open("mlip_input.cfg", "w")

def main():

    lines, energies  = cfg_energy(path)
    n_clusters, labels = cluster(energies)  # Step 1: clustering
    disperse_energy, close_energy = print_clusters(energies, labels)  # Step 2: report + separate clusters
    
    extracted_blocks=extract_data_from_cfg(lines,disperse_energy, close_energy)

    rounded_disperse = [round(float(val), 6) for val in disperse_energy]
    
    close_Energy = []
    
    for j in range(0,len(close_energy)):
    
        if close_energy[j]:
        
           close_Energy.append(random.choice(close_energy[j]))
        
    close_energy_rounded = [round(float(val), 6) for val in close_Energy]
    
    print(len(close_energy_rounded))
    
    n_selected = 0
    n_initial = 0
    for blocks in extracted_blocks:
        
        Energy = blocks["energy"]
        
        Energy_round = round(Energy, 6)
        size = blocks["size"]
        n_initial +=1
        if (Energy_round in rounded_disperse or Energy_round in close_energy_rounded):
           
           n_selected += 1
           
           supercell = blocks["Supercell"]
           atom_infos = blocks["atoms"]
           Stress = blocks["stress"]
           
           write_conf(size,supercell,atom_infos,Stress,Energy)
    
    print('Configurations loaded initially:',n_initial)      
    print('The number of selected configurations:',n_selected)

    #path_to_mlip="/home/karimzongo/mlip-2/bin/mlp mindist /home/karimzongo/Desktop/ETS/SIMULATION/POSTDOC/Na-K-Cl/py/mlip_input.cfg"
    
    #subprocess.call(path_to_mlip, shell=True)
     
if __name__ == "__main__":

    main()




     
     
               
           
             
