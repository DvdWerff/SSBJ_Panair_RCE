import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.ndimage import rotate


def transform_airfoils(airfoils, airfoil_folder_path, root_chord, tip_chord, sweep, half_span, tc_ratio, twist=0, kink_chord=0, sweep_outboard = 0, b_kink = 0):
    '''
    Transform airfoils to have their desired position, chord, and twist.
    Results stored in new csv files with their location names in the "Airfoils" folder.
    :param airfoils: dictionary containing the airfoils used at each location.
                    Currently only supports 'root_airfoil' and 'tip_airfoil'.
    :param root_chord: Chord length of the root airfoil [m]
    :param tip_chord: Chord length of the tip airfoil [m]
    :param sweep: Leading edge sweep angle [deg]
    :param half_span: Half of the total span (Span of a single wing) [m]
    :param tc_ratio: thickness over chord ratio
    :return:
    '''

    # dihedral = 0, twist = 0


    if 'kink_airfoil' in airfoils.keys():
        x_LE_kink = b_kink * np.tan(np.radians(sweep))
        x_LE_tip = x_LE_kink + (half_span - b_kink) * np.tan(np.radians(sweep_outboard))
        twist_kink = b_kink / half_span * twist
    else:
        x_LE_tip = half_span*np.tan(np.radians(sweep))


    for loc, name in airfoils.items():
        if loc == 'root_airfoil':
            transform_airfoil(airfoil_loc=loc,
                              airfoil_folder_path=airfoil_folder_path,
                              airfoil_name=name,
                              x_LE=0,
                              chord=root_chord,
                              tc_ratio=tc_ratio)
        elif loc == 'kink_airfoil':
            transform_airfoil(airfoil_loc=loc,
                              airfoil_folder_path=airfoil_folder_path,
                              airfoil_name=name,
                              x_LE=x_LE_kink,
                              chord=kink_chord,
                              tc_ratio=tc_ratio,
                              twist = twist_kink)
        else:
            transform_airfoil(airfoil_loc=loc,
                              airfoil_folder_path=airfoil_folder_path,
                              airfoil_name=name,
                              x_LE=x_LE_tip,
                              chord=tip_chord,
                              tc_ratio=tc_ratio,
                              twist=twist)

    return
def transform_airfoil(airfoil_loc, airfoil_folder_path, airfoil_name, x_LE, chord, tc_ratio, z_LE = 0, twist= 0):
    """
    Reads an airfoil csv, applies the necessary scaling, translation, rotation (TBD),
        and writes a new csv with the updated coordinates.
    :param name_airfoil: the name of the airfoil
    :param airfoil_loc: spanwise location of the airfoil. For now: either root_airfoil or tip_airfoil
    :param x_LE: x coordinate of the leading edge [m]
    :param chord: chord of the airfoil [m]
    :param z_LE: z coordinate of the leading edge [m]
    :param twist: twist of the airfoil [deg]
    :param tc_ratio: thickness over chord ratio
    :return:
    """

    airfoil_base_csv = os.path.join(airfoil_folder_path, f'{airfoil_name}-formatted.csv')
    baseline_df = pd.read_csv(airfoil_base_csv)
    columns = baseline_df.columns
    coords = baseline_df.values

    #general scaling
    coords_scaled = coords*(chord/100)


    #thickness scaling
    thickness = max(coords_scaled[:, 2] - coords_scaled[:, 4])
    tc_ratio_original = thickness/chord
    coords_scaled[:, [2, 4]] = coords_scaled[:, [2, 4]] * (tc_ratio / tc_ratio_original)

    #rotation
    coords_rotated = rotate_airfoils(coords_scaled, twist)

    # #validate rotation by plotting
    # fig, ax = plt.subplots()
    # ax.scatter(coords_scaled[:, 1], coords_scaled[:, 2], c='b')
    # ax.scatter(coords_scaled[:, 3], coords_scaled[:, 4], c='b')
    #
    # ax.scatter(coords_rotated[:, 1], coords_rotated[:, 2], c='r')
    # ax.scatter(coords_rotated[:, 3], coords_rotated[:, 4], c='r')
    #
    # plt.show()
    #
    # print('keep plot active ')

    #translation
    coords_rotated[:, [1, 3]] = coords_rotated[:, [1, 3]] + x_LE
    coords_rotated[:, [2, 4]] = coords_rotated[:, [2, 4]] + z_LE

    df_transformed_coords = pd.DataFrame(coords_rotated, columns=columns)
    df_transformed_coords.to_csv(os.path.join(airfoil_folder_path, f'{airfoil_loc}.csv'))
    return

def rotate_airfoils(coords, angle):

    def rotate_matrix (arr, angle):
        x = arr[:,0]
        y = arr[:,1]

        angle = np.radians(angle)
        xr = (x * np.cos(angle)) - (y * np.sin(angle))
        yr = (x * np.sin(angle)) + (y * np.cos(angle))

        rot_arr = np.stack((xr,yr))
        return rot_arr.transpose()

    upsurf = coords[:, 1:3]
    lowsurf = coords[:, 3:]

    #negative angle for downwards angle
    upsurf_rot = rotate_matrix(upsurf, -angle)
    lowsurf_rot = rotate_matrix(lowsurf, -angle)

    coords_rot = coords.copy()
    coords_rot[:, 1:3] = upsurf_rot
    coords_rot[:, 3:] = lowsurf_rot

    return coords_rot

if __name__ == '__main__':
    root_airfoil = 'naca64206-il'
    tip_airfoil = 'naca64206-il'
    airfoils = dict(
        root_airfoil=root_airfoil,
        tip_airfoil=tip_airfoil
    )
    airfoil_folder_path = r'C:\Users\Daan\PycharmProjects\SSBJ_Panair_RCE\Airfoils'

    root_chord = 10
    tip_chord = 2.5
    sweep = 55
    half_span = 20
    tc_ratio = 0.05

    twist = 45

    transform_airfoils(airfoils, airfoil_folder_path, root_chord, tip_chord, sweep, half_span, tc_ratio, twist)
