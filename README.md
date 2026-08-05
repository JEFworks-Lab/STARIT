# STARIT

**STARIT (Spatial Transcriptomics As Rasterized Image Tensors)** is an open-source Python package that turns im-SRT data point coordinates into rasterized image tensors. This is the `STARIT` Python documentation
website. Questions, suggestions, or problems should be submitted as
[GitHub issues](https://github.com/JEFworks-Lab/STARIT/issues).

## Overview
<p align="center">

<img src="https://github.com/JEFworks-Lab/STARIT/blob/main/images/0STARIT_fig1.png?raw=true" height="600"/>

</p>

#### Motivation 

Imaging-based spatially resolved transcriptomics (im-SRT) provide high-throughput molecular-resolution spatial characterization of genes within cells. Conventional analysis methods to identify cell-types and states in im-SRT data rely on gene count matrices derived from tallying the number of mRNA molecules detected for each gene per segmented cell, thereby overlooking subcellular heterogeneity that can be useful in defining cell states.  

#### Results 

To take advantage of the molecular-resolution information in im-SRT data and potentially identify cell-states based on subcellular heterogeneity, we developed `STARIT` (Spatial Transcriptomics As Rasterized Image Tensors). `STARIT` converts transcripts within segmented cells in im-SRT data into an image-based tensor representation that can be combined with deep learning computer vision models for downstream analysis. Using simulated data, we demonstrate that `STARIT` distinguishes expected transcriptionally distinct cell types and further separates cell states based on subcellular transcript localization, which conventional gene count analysis fails to capture. Likewise, using real im-SRT data, we demonstrate how `STARIT` can offer comparable results to conventional gene count analysis as well as delineate rotational variation not captured by conventional gene count analysis. By providing a standardized framework to represent subcellular molecular information in im-SRT data, coupled with future technological advancements, `STARIT` will enable deeper insights into subcellular heterogeneity and enhance the identification and characterization of cell types and states that are overlooked by gene count representations. 

## Requirements

This package requires **Python 3.11 or later**.

To ensure compatibility, you must have a version that satisfies the `python_requires` setting in `pyproject.toml`.

## Installation
From GitHub (main branch)
```bash
pip install "git+https://github.com/JEFworks-Lab/STARIT.git#egg=starit"

```

### Test/Install Locally
#### Through `pip`
0) Clean old build junk (if applicable)
```bash
pip install "git+https://github.com/JEFworks-Lab/STARIT.git#egg=starit"

```
1) Fresh virtual environment (venv)
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip

```

2) Editable install
```bash
pip install -e .

```

3) Build the wheel/sdist
```bash
pip install -U build twine
python -m build
twine check dist/*
```

4) Try the built wheel in a *clean* env
```bash
deactivate
python -m venv .venv-clean
source .venv-clean/bin/activate
python -m pip install -U pip
pip install dist/starit-*.whl
```

#### Through Conda

```bash
git clone git@github.com:JEFworks-Lab/STARIT.git
cd STARIT
conda create --name starit python=3.11
pip install -e .
```

Should be able to import and use STARIT like so:

```python 
from starit import starit
```

## Tutorials
- [Running STARIT on Simulated Dataset](https://github.com/JEFworks-Lab/STARIT/blob/main/simulated_data_example.ipynb)
- [Running STARIT on Simulated Dataset with Noise](https://github.com/JEFworks-Lab/STARIT/blob/main/simulated_data_example.ipynb)


## Citation

Our paper describing `STARIT` is availible on **bioRxiv**:

[Velazquez D. et al. (2025), "Spatial Transcriptomics As Rasterized Image Tensors (STARIT) characterizes cell states with subcellular molecular heterogeneity", *bioRxiv*](https://doi.org/10.64898/2025.12.18.695193)

