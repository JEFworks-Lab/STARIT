import numpy as np
import matplotlib.pyplot as plt

def normalize(arr, t_min=0, t_max=1):
    """Linearly normalizes an array between two specifed values.
    
    Parameters
    ----------
    arr : numpy array
        array to be normalized
    t_min : int or float
        Lower bound of normalization range
    t_max : int or float
        Upper bound of normalization range
    
    Returns
    -------
    norm_arr : numpy array
        1D array with normalized arr values
        
    """
    
    diff = t_max - t_min
    diff_arr = np.max(arr) - np.min(arr)
    min_ = np.min(arr)
        
    norm_arr = ((arr - min_)/diff_arr * diff) + t_min
    
    return norm_arr

def get_bounding_box(x, y, expand=1.1):
    ''' Get mininmum and maximum 2D coordinates in all gene molecule positions of your data to create standardized image tensor representations (maintain cell aspecft ratio)
    
    Paramters
    ---------
    x : numpy array of length N
        x location of all molecule points in all cells to serve as segmentation bounds for bounding box
    y : numpy array of length N
        y location of all molecule points in all cells to serve as segmentation bounds for bounding box
    expand : float
        Factor to expand sampled area beyond molecules. Defaults to 1.1.
        
    Returns
    -------
    min_x : float
        Minimum x coordinate of bounding box
    max_x : float
        Maximum x coordinate of bounding box
    min_y : float
        Minimum y coordinate of bounding box
    max_y : float
        Maximum y coordinate of bounding box
        
    '''
    min_x = np.min(x)
    max_x = np.max(x)
    min_y = np.min(y)
    max_y = np.max(y)
    min_x,max_x = (min_x+max_x)/2.0 - (max_x-min_x)/2.0*expand, (min_x+max_x)/2.0 + (max_x-min_x)/2.0*expand
    min_y,max_x = (min_y+max_y)/2.0 - (max_y-min_y)/2.0*expand, (min_y+max_y)/2.0 + (max_y-min_y)/2.0*expand
    
    return min_x, max_x, min_y, max_y

def starit(bounding_box, x, y, g=np.ones(1), dx=1.0, blur=1.0, draw=10000, wavelet_magnitude=False, use_windowing=True):
    ''' Rasterize a spatial transcriptomics dataset into a density image
    
    Paramters
    ---------
    bounding_box: list or tuple of 4 floats
        Bounding box coordinates as [min_x, max_x, min_y, max_y]
    x : numpy array of length N
        x location of molecule points to be rasterized
    y : numpy array of length N
        y location of molecule points to be rasterized
    g : numpy array of length N
        RNA count of cells
        If not given, density image is created
    dx : float
        Pixel size to rasterize data (default 30.0, in same units as x and y)
    blur : float or list of floats
        Standard deviation of Gaussian interpolation kernel.  Units are in 
        number of pixels.  Can be aUse a list to do multi scale.
    draw : int
        If True, draw a figure every draw points return its handle. Defaults to False (0).
    wavelet_magnitude : bool
        If True, take the absolute value of difference between scales for raster images.
        When using this option blur should be sorted from greatest to least.
    
        
    Returns
    -------
    X  : numpy array
        Locations of pixels along the x axis
    Y  : numpy array
        Locations of pixels along the y axis
    M : numpy array
        A rasterized image with len(blur) channels along the first axis
    fig : matplotlib figure handle
        If draw=True, returns a figure handle to the drawn figure.
        
    Raises
    ------    
    Exception 
        If wavelet_magnitude is set to true but blur is not sorted from greatest to least.
    '''
    
    # set blur to a list
    if not isinstance(blur, list):
        blur = [blur]
    nb = len(blur)
    blur = np.array(blur)
    n = len(x)
    maxblur = np.max(blur)  # for windowing
    
    if wavelet_magnitude and np.any(blur != np.sort(blur)[::-1]):
        raise Exception('When using wavelet magnitude, blurs must be sorted from greatest to least')
    
    # Unpack bounding box coordinates
    min_x, max_x, min_y, max_y = bounding_box
    
    X_ = np.arange(min_x, max_x, dx)
    Y_ = np.arange(min_y, max_y, dx)
    
    X = np.stack(np.meshgrid(X_, Y_))  # note this is xy order, not row col order

    W = np.zeros((X.shape[1], X.shape[2], nb))

    if draw: 
        fig, ax = plt.subplots()
    count = 0
    
    g = np.resize(g, x.size)
    if not (g == 1.0).all():
        g = normalize(g)
    
    for x_, y_, g_ in zip(x, y, g):
        if not use_windowing:  # legacy version
            k = np.exp(- ((X[0][..., None] - x_)**2 + (X[1][..., None] - y_)**2) / (2.0*(dx*blur*2)**2))
            k /= np.sum(k, axis=(0, 1), keepdims=True)
            k *= g_
            if wavelet_magnitude:
                for i in reversed(range(nb)):
                    if i == 0:
                        continue
                    k[..., i] = k[..., i] - k[..., i-1]
            W += k
        else:  # use a small window
            r = int(np.ceil(maxblur*4))
            col = np.round((x_ - X_[0])/dx).astype(int)
            row = np.round((y_ - Y_[0])/dx).astype(int)
            
            row0 = np.floor(row-r).astype(int)
            row1 = np.ceil(row+r).astype(int)                    
            col0 = np.floor(col-r).astype(int)
            col1 = np.ceil(col+r).astype(int)
            # we need boundary conditions
            row0 = np.minimum(np.maximum(row0, 0), W.shape[0]-1)
            row1 = np.minimum(np.maximum(row1, 0), W.shape[0]-1)
            col0 = np.minimum(np.maximum(col0, 0), W.shape[1]-1)
            col1 = np.minimum(np.maximum(col1, 0), W.shape[1]-1)
            
            k = np.exp(- ((X[0][row0:row1+1, col0:col1+1, None] - x_)**2 + (X[1][row0:row1+1, col0:col1+1, None] - y_)**2) / (2.0*(dx*blur*2)**2))
            k /= np.sum(k, axis=(0, 1), keepdims=True)  
            k *= g_
            if wavelet_magnitude:
                for i in reversed(range(nb)):
                    if i == 0:
                        continue
                    k[..., i] = k[..., i] - k[..., i-1]
            W[row0:row1+1, col0:col1+1, :] += k

        if draw:
            if not count % draw or count == (x.shape[0]-1):
                print(f'{count} of {x.shape[0]}')
                ax.cla()
                toshow = W - np.min(W, axis=(0, 1), keepdims=True)
                toshow = toshow / np.max(toshow, axis=(0, 1), keepdims=True)
                
                if nb >= 3:
                    toshow = toshow[..., :3]
                elif nb == 2:
                    toshow = toshow[..., [0, 1, 0]]
                elif nb == 1:
                    toshow = toshow[..., [0, 0, 0]]
                
                ax.imshow(np.abs(toshow))
                ax.invert_yaxis()
                ax.axis('off')
                fig.canvas.draw()

        count += 1
        
    W = np.abs(W)
    # we will permute so channels are on first axis
    W = W.transpose((-1, 0, 1))
    extent = (X_[0], X_[-1], Y_[0], Y_[-1])
    
    # rename
    X = X_
    Y = Y_
    if draw:
        output = X, Y, W, fig
    else:
        output = X, Y, W
    return output



def starit_gene(bounding_box, x_gene_pos, y_gene_pos, dx=1, blur=1):
    """Rasterize gene-specific molecular positions of a cell and returns gene rasterized image figure
    
    Paramters
    ---------
    bounding_box: list or tuple of 4 floats
        Bounding box coordinates as [min_x, max_x, min_y, max_y]
    x_gene_pos : numpy array of length N
        x location of coordinate points of gene of interest to be rasterized
    y_gene_pos: numpy array of length N
        y location of coordinate points of gene of interest to be rasterized
    dx : float
        Pixel size to rasterize data (default 30.0, in same units as x and y)
    blur : float or list of floats
        Standard deviation of Gaussian interpolation kernel.  Units are in 
        number of pixels.  Can be aUse a list to do multi scale.
    
    Returns
    -------
    gene_img : matplotlib figure handle
        Returns a figure handle to the drawn figure.
        
    """
    _, _, _, gene_img = starit(bounding_box, x_gene_pos, y_gene_pos, dx=dx, blur=blur)
    return gene_img