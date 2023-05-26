import pandas as pd
import numpy as np

def format_csv(filename):
    """
    Function to format airfoil coordinates downloaded from airfoiltools.com to the
    :param filename: filename (without extension) of the csv file to be formatted
    :return: creates formatted csv in the same directory
    """

    with open(f"{filename}.csv", 'r') as f:
        df = pd.read_csv(f)

    vals = df.values
    coords = vals[8:59]

    len_per_surf = len(coords)/2

    surf_up_flipped = coords[:int(np.ceil(len_per_surf))].astype(float)
    surf_up = np.flip(surf_up_flipped,0)
    surf_low = coords[int(np.floor(len_per_surf)):].astype(float)

    formatted_array = np.zeros((int(np.ceil(len_per_surf)),4))
    formatted_array[:,:2] = surf_up
    formatted_array[:,2:] = surf_low

    columns = ["xup", "zup", "xlow", "zlow"]

    formatted_df = pd.DataFrame(formatted_array, columns=columns)
    formatted_df.to_csv(f"{filename}-formatted.csv")
    return

if __name__ == '__main__':
    airfoil = 'naca64206-il'
    format_csv(airfoil)

    #validate
    print(pd.read_csv(f'{airfoil}-formatted.csv'))


