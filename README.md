# cellLandscape

## Overview

<p>
  <img src="overview.png" width=50% height=50% />
</p>

## Installation
Dependencies 
* Python >= 3.6, anndata 0.7.6, numpy 1.22.4, scipy 1.7.1, tqdm 4.64.0, scanpy 1.8.1

You can clone the git repository by, 

```
git clone https://github.com/CompCy-lab/cellLandscape.git
```

## Example usage
To compute the kernel mean embedding, first read in a preprocessed `.h5ad` object. This dataset contains multiple profiled single-cell samples. 

```python
import scanpy as sc
adata = sc.read_h5ad('nk_cell_preprocessed.h5ad')
```
Then simply compute the kernel mean embedding as,

```python
# Inputs:
# adata: annotated data object (dimensions = cells x features)
# sample_set_key: string referring to the key within adata.obs that contains the samples to compute the embedding
#   ~ If sample_set_key is None, will use all cells as a single sample 
# gamma: scale for standard deviation of the normal distribution within random Fourier frequency feature computation  
# D: dimensionality of the random Fourier frequency features, D/2 sin and D/2 cos basis 
# frequency_seed: random state parameter 
# -----------------------
    
# Returns:
# X_embed: mean embedding (dimensions = samples x D)

# -----------------------
from embed import mean_embed
X_embed = mean_embed(adata, sample_set_key = 'FCS_File', gamma = 1, D = 2000, frequency_seed = 0)
```
