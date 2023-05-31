# This folder has workflow templates for spatial transcriptomics data.

* spatial_10X.ipynb: workflow for spatial cell deconvolution.
* cell2location: workflow for cell type deconvolution using cell2location.

In folder projects, it has jupyter notebooks for each separate customized analysis.

## Docker container download in HPC

For workflow in spatial_10X.ipynb, please use the following container.

	singularity pull docker://shl198/sc_ppl:202306

For workflow in cell2location.ipynb, please use the following container, which can be directly download:
	
	wget https://cell2location.cog.sanger.ac.uk/singularity/cell2location-v0.06-alpha.sif

This webpage has the tutorial to set up environment to run cell2location. https://cell2location.readthedocs.io/en/latest/dockersingularity.html.

## How to run the workflow
Users can follow the steps in file run_jupyterlab_in_HPC.md in parent directory to setup the environment. The detailed steps are described in the jupyter notebooks. 
