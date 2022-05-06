from backspinpy import CEF_obj
import anndata
import numpy as np
import scanpy as sc
import pandas as pd


cef = CEF_obj()
cef.readCEF("LizardA_expression.cef")
df = pd.DataFrame(data = cef.matrix, index=cef.row_attr_values[0], columns=cef.col_attr_values[0])
df = df.transpose()
#print(df.head(5))
count= df.iloc[:,:]
adata = anndata.AnnData(count)
adata.obs_names = df.index
adata.var_names = df.columns
print(adata.var_names)

sc.pl.highest_expr_genes(adata, n_top=20, )

#filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=7)
#violin plot
sc.pl.violin(adata, ['n_genes'],
             jitter=0.4, multi_panel=True)

sc.pp.normalize_total(adata, target_sum=1e4) #normalize every sample to 10000 counts
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

#PCA
sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca(adata, color='LAMP5')
sc.pl.pca_variance_ratio(adata, log=True)

print(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=18)

#sc.tl.paga(adata)
#sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata)


#sc.tl.leiden(adata, key_added = "leiden_1.0") # default resolution in 1.0
#sc.tl.leiden(adata, resolution = 0.6, key_added = "leiden_0.6")
#sc.tl.leiden(adata, resolution = 0.4, key_added = "leiden_0.4")
sc.tl.leiden(adata, resolution = 1, key_added = "leiden")
#sc.pl.umap(adata,color = 'leiden_1.4')
#sc.pl.umap(adata, color=['leiden_0.4', 'leiden_0.6', 'leiden_1.0','leiden_1.4'])

sc.tl.louvain(adata, resolution = 1,key_added = "louvain")
sc.pl.umap(adata,color = ['louvain','leiden','RTCD1','PVALB'])
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon') #t-test, logistic regression (log-reg)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
print(adata.uns['rank_genes_groups']['names']['0'])