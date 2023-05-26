import numpy as np
from scipy.ndimage import rotate
from matplotlib import pyplot as plt

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
    arr = np.array([[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00],
                    [1.0000e+00, 4.5900e-01, 5.4200e-01, 5.4100e-01, -4.4200e-01],
                    [2.0000e+00, 7.0400e-01, 6.6400e-01, 7.9600e-01, -5.2400e-01],
                    [3.0000e+00, 1.1980e+00, 8.5900e-01, 1.3020e+00, -6.4500e-01]], dtype=float)

    arr_rot = rotate_airfoils(arr, 10)

    fig, ax = plt.subplots()
    ax.scatter(arr[:, 1], arr[:, 2], c='b')
    ax.scatter(arr[:, 3], arr[:, 4], c='b')

    ax.scatter(arr_rot[:, 1], arr_rot[:, 2], c='r')
    ax.scatter(arr_rot[:, 3], arr_rot[:, 4], c='r')

    plt.show()

    print('keep plot active ')

