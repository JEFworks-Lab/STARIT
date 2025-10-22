import numpy as np
from starit.starit import starit

def test_starit_shapes():
    x_cell_seg = np.array([0, 100])
    y_cell_seg = np.array([0, 50])
    rng = np.random.default_rng(0)
    x = rng.uniform(10, 90, 100)
    y = rng.uniform(10, 40, 100)

    X, Y, M, fig = starit(x_cell_seg, y_cell_seg, x, y, dx=2.0, blur=1.0, expand=1.1, draw=10000)
    assert M.ndim == 3  # (C, H, W)
    assert M.shape[0] == 1
    assert len(X) > 0 and len(Y) > 0
    assert M.shape[1] == len(Y) and M.shape[2] == len(X)