# STARIT

**STARIT (Spatial Transcriptomics As Rasterized Image Tensors)** turns im-SRT data point coordinates into rasterized image tensors.

## Requirements

This package requires **Python 3.11 or later**.

To ensure compatibility, you must have a version that satisfies the `python_requires` setting in `pyproject.toml`.

## Install
Eventually:
```bash
pip install starit
```

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

