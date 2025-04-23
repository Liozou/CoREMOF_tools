import pandas as pd

node_df = pd.read_csv('MOFid-v5_node.csv', header=0)

node_dict = {}

with open('MOFid-v5_node.csv', 'r') as f:
    for index, line in enumerate(f):
        if index >= 1:
            data = line.strip().split(',')
            mofname = data[1] + '_' + data[2] + '_pacman'
            nodename = f'[{data[3]}_{data[4]}]'
            if mofname not in node_dict:
                node_dict[mofname] = []
                node_dict[mofname].append(nodename)
            else:
                node_dict[mofname].append(nodename)

from ase.io import read as ase_read
from ase.io import write as ase_write
import networkx as nx
import selfies as sf

def dict2str(dct):
    """Convert symbol-to-number dict to str.
    """
    return ''.join(symb + (str(n)) for symb, n in dct.items())

metals  = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
          'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi', 'Po',
          'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
          'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
          'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U', 'Tm', 'Yb', 'Lu',
          'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
         ]


test_mofid = '[Ni].[O]P(=O)(Cc1ccc(cc1)C1=NC(=[N]=C([N]1)c1ccc(cc1)CP(=O)(O)O)c1ccc(cc1)CP(=O)(O)[O])O.[O]P(=O)(Cc1ccc(cc1)C1=NC(=[N]=C([N]1)c1ccc(cc1)CP(=O)(O)O)c1ccc(cc1)CP(=O)([O])O)O MOFid-v1.ERROR.cat0;ABAVIM01_ASR'
def convert_mofid_to_mofidv2(mofid):
    def join_mofidv2(selfies, node_type, v1_details):
        v2_details = v1_details.replace("v1", "v2")
        
        if selfies == '*' and node_type == '*':
            return '*' + ' ' + v2_details
        else:
            return node_type + '.' + selfies + ' ' + v2_details
    
    def contain_metal(smile):
        check = False
        for m in metals:
            if m in smile:
                check = True
        return check
    
    def smiles_to_selfies(smiles):    
        selfies = []
        for smile in smiles:
            if contain_metal(smile):
                continue
            else:
                try:
                    selfie = sf.encoder(smile)
                except:
                    selfie = 'ERROR'
                selfies.append(selfie)
        return '.'.join(selfies)

    smiles = mofid.split()[0].split('.')
    v1_details = mofid.split()[1]

    if smiles == ['*']:
        selfies = '*'
        node_type = '*'
    else:
        selfies = smiles_to_selfies(smiles)
        mofname = v1_details.split(';')[-1]
        try:
            nodes = node_dict[mofname]
        except KeyError:
            print(mofname)
            nodes = ['*']
        if len(nodes) == 1:
            node_type = nodes[0]
        else:
            node_type = '.'.join(nodes)
        
    mofid_v2 = join_mofidv2(selfies, node_type, v1_details)
    return mofid_v2
convert_mofid_to_mofidv2(test_mofid) 



# ASR
mofids = []
mofids_v2 = []
with open('ASR_mofid-v1.txt', 'r') as f:
    for line in f:
        mofids.append(line.strip())
with open('ASR_mofid-v2.txt', 'w') as f:        
    for idx, mofid in enumerate(mofids):
        mofid_v2 = convert_mofid_to_mofidv2(mofid)
        f.write('{}\n'.format(mofid_v2))



# FSR
mofids = []
mofids_v2 = []
with open('FSR_mofid-v1.txt', 'r') as f:
    for line in f:
        mofids.append(line.strip())
with open('FSR_mofid-v2.txt', 'w') as f:        
    for idx, mofid in enumerate(mofids):
        mofid_v2 = convert_mofid_to_mofidv2(mofid)
        f.write('{}\n'.format(mofid_v2))


# Ions
mofids = []
mofids_v2 = []
with open('ion_mofid-v1.txt', 'r') as f:
    for line in f:
        mofids.append(line.strip())
with open('ion_mofid-v2.txt', 'w') as f:        
    for idx, mofid in enumerate(mofids):
        mofid_v2 = convert_mofid_to_mofidv2(mofid)
        f.write('{}\n'.format(mofid_v2))



refs1 = []
refs2 = []
linkers = []
topos = []
cats = []
metal_nodes = []
metals  = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
          'Al', 'Ga', 'Ge', 'In', 'Sn', 'Sb', 'Tl', 'Pb', 'Bi', 'Po',
          'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
          'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
          'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
          'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U', 'Tm', 'Yb', 'Lu',
          'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
         ]

def contain_metal(smile):
    check = False
    for m in metals:
        if m in smile:
            check = True
    return check

with open('all_mofid-v1.txt', 'r') as f:
    for line in f:
        if line.startswith("*"):
            continue
        else:
            data = line.strip().split(';')
            refs1.append(data[1])
            v1 = data[0].split()
            big_smile = v1[0].split('.')
            others = v1[1]
            
            final_smile = []
            others = others.split('.')

            for smile in big_smile:
                if contain_metal(smile):
                    continue
                else:
                    final_smile.append(smile)s
            linker = '.'.join(final_smile) #Linker SMILES
            topo = others[1]
            cat = others[-1]
            topos.append(topo)
            cats.append(cat)
            linkers.append(linker)
            
with open('all_mofid-v2.txt', 'r') as f:
    for line in f:
        if line.startswith("*"):
            continue
        else:
            data = line.strip().split(';')
            refs2.append(data[1])
            v2 = data[0].split()
            big_smile = v2[0].split('.')
            others = v2[1]
            
            final_smile = []

            for smile in big_smile:
                if contain_metal(smile):
                    final_smile.append(smile)
            metal_node = '.'.join(final_smile) # Metal node
            print(metal_node)
            metal_nodes.append(metal_node)
with open('data.csv', 'w') as f:
    for i, j, k, l in zip(refs1, linkers, topos, cats):
        try:
            id_ref2 = refs2.index(i)
            #print(id_ref2)
            metal_node = metal_nodes[id_ref2]
        except ValueError:
            metal_node = None
        f.write(f"{i}|{metal_node}|{j}|{k}|{l}\n")