import matplotlib.pyplot as plt

def save_image(fig, filename, format='png'):
    file = filename
    fig.savefig(file, bbox_inches='tight', pad_inches=0, format=format)