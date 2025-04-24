"""
[1] Men et al., 2021. https://doi.org/10.1016/j.solener.2021.07.051.
"""
from niopy.geometric_transforms import tgs2tube, ang
from numpy import array, sqrt, rad2deg
from pandas import read_csv, DataFrame


########################################################################################################################
# Functions for this program ###########################################################################################

def theoretical_acceptance(mirror_width: float, center_distance: float, n_mirrors: int,
                           receiver_height: float, tube_radius: float):

    total_width = center_distance * (n_mirrors - 1) + mirror_width

    f2 = array([total_width / 2, 0])

    tube_center = array([0, receiver_height])

    t1, t2 = tgs2tube(point=f2, tube_radius=tube_radius, tube_center=tube_center)
    t1 = t1 if t1[0] < t2[0] else t2

    half_acceptance = rad2deg(ang(t1 - f2, array([0, 1]))).round(2)

    return half_acceptance


########################################################################################################################
########################################################################################################################


nbr_mirrors = 25

# importing raw data
Men2021_data = read_csv('Men2021_raw_data.csv')

# Analysis of CPC truncation ###########################################################################################
"""
According to the model presented by the authors, for a full unrolled CPC (no truncation), the sum of the acceptance
half-angle (theta_a) and the max unrolled angle (theta_{max}) should be 270º.

Then it creates a new column where true indicates that the optimum cpc geometry reported by the authors is full unrolled 
(TRUE) or not (FALSE).
"""

Men2021_data['check_full_cpc'] = [(Men2021_data['theta_a'][i] + Men2021_data['theta_max'][i]) == 270.0
                                  for i in range(Men2021_data['theta_a'].shape[0])]
########################################################################################################################

# Analysis of cpc acceptance angle #####################################################################################

Men2021_data['theoretical_acceptance'] = [theoretical_acceptance(mirror_width=w,
                                                                 center_distance=d,
                                                                 n_mirrors=nbr_mirrors,
                                                                 receiver_height=h,
                                                                 tube_radius=ra)
                                          for w, d, h, ra in zip(Men2021_data['W'],
                                                                 Men2021_data['D'],
                                                                 Men2021_data['H_R'],
                                                                 Men2021_data['r_a'])
                                          ]

Men2021_data['acceptance_diff'] = Men2021_data['theta_a'] - Men2021_data['theoretical_acceptance']

########################################################################################################################

# Analysis of Iparraguirre design ######################################################################################
# n_steps = (nbr_mirrors - 1) / 2 if nbr_mirrors % 2 != 0 else nbr_mirrors/2
#
# Men2021_data['Iparraguirre design'] = [sqrt((n_steps * d)**2 + h**2)
#                                        if s == 'Parabolic'
#                                        else 2 * sqrt((n_steps * d)**2 + h**2)
#                                        for d, h, s in zip(Men2021_data['D'],
#                                                           Men2021_data['H_R'],
#                                                           Men2021_data['mirror_shape'])
#                                        ]
#
# Men2021_data['Iparraguirre design'] = Men2021_data['Iparraguirre design'].round(2)
#
#
# Men2021_data['curvature_diff'] = [abs(c2 - c1)
#                                   if s != 'Flat'
#                                   else 'NA'
#                                   for c1, c2, s in zip(Men2021_data['MSF_i'],
#                                                        Men2021_data['Iparraguirre design'],
#                                                        Men2021_data['mirror_shape'])
#                                   ]

# Exporting a new data frame with columns to indicate the new analysis done
Men2021_data.to_csv('Men2021_analyzed_data.csv', index=False)

########################################################################################################################
########################################################################################################################
