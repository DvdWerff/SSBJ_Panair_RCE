import os
import numpy as np
import pandas as pd

"""
Script to convert .dat files to the .csv format required for pyPanair
"""
files = os.listdir(os.getcwd())

for file in files:
    if '.dat' not in file:
        continue
    with open(file, 'r') as f:
        lines = f.read()

    filename = file.replace('.dat','')

    lines = lines.split('\n')
    coords = np.zeros((len(lines),2))
    for i, line in enumerate(lines):
        coords[i] = [float(val) for val in line.split()]

    df = pd.DataFrame(coords)
    df.to_csv(f'{filename}.csv', index=False)
