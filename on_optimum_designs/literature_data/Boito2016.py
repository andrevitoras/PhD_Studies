"""
Boito, P., Grena, R., 2016. Optimization of the geometry of Fresnel linear collectors.
Solar Energy 135, 479–486. https://doi.org/10.1016/j.solener.2016.05.060
"""
from numpy import array

if __name__ == '__main__':

    # Curved primary field with a flat primary field.

    # Base value used for the optimization that are kept constant (Section 5.1, p. 483).
    number_primaries = 25
    receiver_height = 10.0
    receiver_width = 0.4
    ####################################################################################################################
    ####################################################################################################################

    # For the UC #######################################################################################################
    # Optimum configuration considering uniform width, shift, and curvature.

    UC_s = 1.21
    UC_w = 2 * 0.5
    UC_R = 2 * 16.85

    s_to_w = UC_s / UC_w
    Wp = (number_primaries - 1) * UC_s + UC_w
    Hr_to_Wp = receiver_height / Wp
    FF = number_primaries * UC_w / Wp
    H_to_nw = receiver_height / (number_primaries*UC_w)

    print(f'Boito UC configuration has: '
          f's/w = {round(s_to_w, 3)}, '
          f'Hr/wp = {round(Hr_to_Wp, 3)}, '
          f'FF = {round(FF, 3)}, '
          f'H/nw = {round(H_to_nw, 3)}')

    ####################################################################################################################
    ####################################################################################################################

    # For OP case ######################################################################################################
    # Optimum configuration considering uniform width, uniform curvature, and variable shift.

    OP_w = 2 * 0.5
    OP_R = 2 * 16.87
    OP_x = [0., 1.14, 2.29, 3.43, 4.57, 5.74, 6.87, 8.03, 9.27, 10.51, 11.73, 13.07, 14.46]

    Wp = 2 * ((OP_x[-1] - OP_x[0]) + 0.5 * OP_w)
    Hr_to_Wp = receiver_height / Wp
    FF = number_primaries * OP_w / Wp
    H_to_nw = receiver_height / (number_primaries*OP_w)

    print(f'Boito OP configuration has: '
          f's/w = var, '
          f'Hr/wp = {round(Hr_to_Wp, 3)}, '
          f'FF = {round(FF, 3)}, '
          f'H/nw = {round(H_to_nw, 3)}')

    ####################################################################################################################
    ####################################################################################################################

    # For OF case ######################################################################################################
    # Optimum configuration considering: uniform width, uniform shift, and variable curvature.

    OF_s = 1.29
    OF_w = 2 * 0.53
    OF_R = 2 * array([13.13, 12.92, 13.07, 13.38, 13.70, 14.04, 15.04, 16.33, 17.50, 18.68, 19.94, 21.18, 22.41])

    s_to_w = OF_s/OF_w
    Wp = (number_primaries - 1) * OF_s + OF_w
    Hr_to_Wp = receiver_height / Wp
    FF = number_primaries * OF_w / Wp
    H_to_nw = receiver_height / (number_primaries*OF_w)

    print(f'Boito OF configuration has: '
          f's/w = {round(s_to_w, 3)}, '
          f'Hr/wp = {round(Hr_to_Wp, 3)}, '
          f'FF = {round(FF, 3)}, '
          f'H/nw = {round(H_to_nw, 3)}')

    ####################################################################################################################
    ####################################################################################################################

    # For NU case ######################################################################################################
    # Optimum configuration considering: variable width, shift, and curvature.

    NU_x = [0., 1.51, 3.97, 6.3, 8.02, 9.41, 10.61, 11.65, 12.61, 13.51, 14.35, 15.16, 16.15]
    NU_w = 2 * array([0.37, 0.88, 1.24, 0.82, 0.62, 0.51, 0.44, 0.37, 0.35, 0.32, 0.27, 0.25, 0.44])
    NU_R = 2 * array([15.97, 10.56, 12.05, 13.83, 15.41, 16.78, 18.07, 19.04, 20.05, 21.02, 21.93, 22.49, 22.71])

    Wp = 2 * ((NU_x[-1] - NU_x[0]) + 0.5 * NU_w[-1])
    Hr_to_Wp = receiver_height / Wp
    FF = 2*NU_w.sum() / Wp
    H_to_nw = receiver_height / (2*NU_w.sum())

    print(f'Boito NU configuration has: '
          f's/w = var, '
          f'Hr/wp = {round(Hr_to_Wp, 3)}, '
          f'FF = {round(FF, 3)}, '
          f'H/nw = {round(H_to_nw, 3)}')

    ####################################################################################################################
    ####################################################################################################################
