from mofid.run_mofid import cif2mofid

def run_v1(structure):
    mofid_v1 = cif2mofid(structure)['mofid']
    return mofid_v1