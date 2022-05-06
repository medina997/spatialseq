import scanpy as sc
import numpy as np
import h5py
filename = "Lizard_gene_values.hdf5"

hf = h5py.File("Lizard_gene_values.hdf5", 'r')
print(hf.keys())
n1 = hf.get('genes')
print(n1)
n1 = np.array(n1)
print(n1)

n2 = hf.get('coordinates')
print(n2)
n2 = np.array(n2)
print(n2)


""""
with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    # Get the data
    data = list(f[a_group_key])
adata = sc.read_hdf("Lizard_gene_values.hdf5",'coordinates')
"""
