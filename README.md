# STARIT

**STARIT (Spatial Transcriptomics As Rasterized Image Tensors)** turns im-SRT data point coordinates into rasterized image tensors.

## Install
Eventually:
```bash
pip install starit
```

From GitHub (main branch)
```bash
pip install "git+https://github.com/JEFworks-Lab/STARIT.git#egg=starit"

```

### Test/install locally
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

Should be able to import and use STARIT like so:

```python 
from starit import starit
```

