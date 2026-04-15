import numpy as np


def a2r(eulerangles):
    ez = eulerangles[0] #kappa
    ey = eulerangles[1] #phi
    ex = eulerangles[2] #omega

    c1 = np.cos(ey)
    c2 = np.cos(ex)
    c3 = np.cos(ez)
    s1 = np.sin(ey)
    s2 = np.sin(ex)
    s3 = np.sin(ez)

    r = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    # r[0] = c1 * c3 - s1 * s2 * s3
    # r[1] = c2 * s3
    # r[2] = s1 * c3 + c1 * s2 * s3

    # r[3] = -c1 * s3 - s1 * s2 * c3
    # r[4] = c2 * c3
    # r[5] = -s1 * s3 + c1 * s2 * c3

    # r[6] = -s1 * c2
    # r[7] = -s2
    # r[8] = c1 * c2

    r[0] = c1 * c3
    r[1] = c1 * s3
    r[2] = -s1
    r[3] = s2 * s1 * c3 - c2 * s3
    r[4] = s2 * s1 * s3 + c2 * c3
    r[5] = s2 * c1
    r[6] = c2 * s1 * c3 + s2 * s3
    r[7] = c2 * s1 * s3 - s2 * c3
    r[8] = c2 * c1

    rot = np.array([[r[0], r[1], r[2]], [r[3], r[4], r[5]], [r[6], r[7], r[8]]])
    return rot