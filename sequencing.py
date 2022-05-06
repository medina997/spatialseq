import pandas as pd
import anndata
import numpy as np
import scanpy as sc
import louvain

"""
df = pd.read_csv("Mouse9_expression.csv")
print(df.head(5))
df1 = df.iloc[:,1:]
counts = df1.to_numpy()
adata = anndata.AnnData(counts)
print(adata)
adata.obs_names = df1.index
adata.var_names = df1.columns
"""
adata = sc.read_10x_mtx(
    'filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
print(adata.obs_names)
print(adata.var_names)
print(adata.uns)

sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
#adata.uns["neighbors"] = [1, 2, 3]

#sc.tl.paga(adata)
#sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata)

#sc.tl.leiden(adata, key_added = "leiden_1.0") # default resolution in 1.0
#sc.tl.leiden(adata, resolution = 0.6, key_added = "leiden_0.6")
#sc.tl.leiden(adata, resolution = 0.4, key_added = "leiden_0.4")
sc.tl.leiden(adata, resolution = 1.4, key_added = "leiden_1.4")
sc.pl.umap(adata,color = 'leiden_1.4')
#sc.pl.umap(adata, color=['leiden_0.4', 'leiden_0.6', 'leiden_1.0','leiden_1.4'])

sc.tl.louvain(adata, resolution = 1.4,key_added = "louvain_1.4")
sc.pl.umap(adata,color = 'louvain_1.4')


