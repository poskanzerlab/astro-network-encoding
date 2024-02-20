# astro-network-encoding - Python scripts

Accompanying the paper "Network-level encoding of local neurotransmitters in cortical astrocytes" (DOI: TODO)

## Prerequisites

Code written and run using Python 3.8.18 (Max Collard).

See `requirements.txt` for versions of other dependencies used.

## Instructions

1. Install `conda`; see instructions [here](https://docs.anaconda.com/free/miniconda/miniconda-install/).
2. Create a `conda` environment using the specification generated from the analysis environment used for the paper, provided in `environments/astro-network-encoding.yml`, for example:
```
conda env create -n astro-network-encoding -f=environments/astro-network-encoding.yml
```
3. Use your prefered method of running Jupyter notebooks (e.g., [VS Code](https://code.visualstudio.com/docs/datascience/jupyter-notebooks)) to run the provided notebooks for the desired figure panel inside of this environment.

### Event tables (`data`)

These scripts take as input HDF-5 files (or `mat`-files in a HDF-compatible format) that contain summarized tables of AQuA event characteristics, in order to  save memory and hard disk space.

The event files used for generating figures in the paper are located in the [Dryad repository](TODO), as `TODO`. The configuration files in the `config` folder are set up to use these files by default. To use this original data:
1. Place the extracted contents of `TODO` into `data/events` (so that `data/events` should contain many subdirectories, such as `X`, `Y`, etc.).

To generate your own event tables from full AQuA `res` files, run `aqua_export.m` on your dataset (see documentation there).

### Ramping cells (`intermediates`)

These scripts also take as input files specifying the "ramping" characteristics of individual cells (see paper *Methods*).

The ramping cell specification files used in the paper are located in the [Dryad repository](TODO), as `TODO`. The configuration files in the `config` folder are set up to use these files by default. To use this original data:
1. Place the extracted contents of `TODO` into `intermediates/ramping` (so that TODO).

To generate your own ramping cell tables from extracted event tables, see TODO.

### Notebooks for figure panels

| Figure Panel | Notebook |
| ------------ | -------- |
| Fig. 3j-l | `net_astro-3j_3k_3l.ipynb` |
| Fig. 4b | `net_astro-4b_s4a_s4b.ipynb` |
| Fig. 4d | `net_astro-4d.ipynb` |
| Fig. 4f | `net_astro-4f_4j_s6b_s6c.ipynb` |
| Fig. 4h | `net_astro-4h_4k_4l_s6e_s6f.ipynb` |
| Fig. 4i | `net_astro-4i.ipynb` |
| Fig. 4j | `net_astro-4f_4j_s6b_s6c.ipynb` |
| Fig. 4k | `net_astro-4h_4k_4l_s6e_s6f.ipynb` |
| Fig. 4l | `net_astro-4h_4k_4l_s6e_s6f.ipynb` |
| Fig. 4m | `net_astro-4m_s5k.ipynb` |
| Ext. Data Fig. 3j | `net_astro-s3j_s5c_s5d_s5e.ipynb` |
| Ext. Data Fig. 4a | `net_astro-4b_s4a_s4b.ipynb` |
| Ext. Data Fig. 4b | `net_astro-4b_s4a_s4b.ipynb` |
| Ext. Data Fig. 5c | `net_astro-s5b.ipynb` |
| Ext. Data Fig. 5c | `net_astro-s3j_s5c_s5d_s5e.ipynb` |
| Ext. Data Fig. 5d | `net_astro-s3j_s5c_s5d_s5e.ipynb` |
| Ext. Data Fig. 5e | `net_astro-s3j_s5c_s5d_s5e.ipynb` |
| Ext. Data Fig. 5k | `net_astro-4m_s5k.ipynb` |
| Ext. Data Fig. 6b | `net_astro-4f_4j_s6b_s6c.ipynb` |
| Ext. Data Fig. 6c | `net_astro-4f_4j_s6b_s6c.ipynb` |
| Ext. Data Fig. 6e | `net_astro-4h_4k_4l_s6e_s6f.ipynb` |
| Ext. Data Fig. 6f | `net_astro-4h_4k_4l_s6e_s6f.ipynb` |
| Ext. Data Fig. 7c | `net_astro-s7c.ipynb` |

## Notes

* Docstrings within individual notebooks are incomplete (marked with "TODO") as of initial release-time; check back here for updates that fill them in!
* There may be subtle differences in p-values or confidence intervals (particularly ones involving resampling) from the exact values for analyses reported in the paper, due to differences in the cached random seeds.
* Additionally, certain outputs seen in notebooks were subsequently produced with lower numbers of iterations than the runs used in the original paper during preparation for distribution, as for certain analyses, full runs with the original number of iterations takes a substantial amount of time; however, numbers of resampling iterations are provided in the corresponding `config` files for users to manipulate.