import matplotlib.pyplot as plt
import seaborn
from matplotlib.ticker import EngFormatter
from numpy import arange, ones, pi, cos, where
from pandas import read_csv, DataFrame, concat
from scipy.interpolate import interp1d

seaborn.set(style='whitegrid')
seaborn.set_context('notebook')

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'


def interpolated_iam_function(file: str):

    k_data = read_csv(file, names=['theta', 'k'])
    k_data['k'] = where(k_data['k'] > 0, k_data['k'], 0)
    fixed_points = DataFrame({'theta': [0., 90.],
                              'k': [1., 0.]
                              })

    df = concat([k_data, fixed_points], join='inner').sort_values(by='theta')

    k_f = interp1d(*df.values.T, kind='cubic')

    return k_f

# PTC Incidence Angle Modifiers #################################
# ptc_kt = ones(shape=angles.shape[0])
#
# ptc_kl_data = read_csv('LFC_transversal.csv', names=['theta', 'k'])
# ptc_kl_data['k'] = where(ptc_kl_data['k'] > 0, ptc_kl_data['k'], 0)
# df = DataFrame({'theta': [0., 90.],
#                 'k': [1., 0.]
#                 })
#
# a = concat([ptc_kl_data, df], join='inner').sort_values(by='theta')
#
# ptc_kl_inter_function = interp1d(*a.values.T, kind='cubic')


def ptc_kl_f(x):
    value = cos(x * pi / 180) - 0.000525 * x - 0.0000286 * x ** 2
    kl = where(value > 0, value, 0)

    return kl


#####################################################################

# LFC Incidence Angle Modifiers #################################

lfc_kt_f = interpolated_iam_function('LFC_transversal.csv')
lfc_kl_f = interpolated_iam_function('LFC_longitudinal.csv')


#####################################################################

if __name__ == '__main__':
    fig_width = 8
    fig_height = 6
    font_size = 8

    lw = 0.75
    ms = 2

    angles = arange(start=0, stop=95, step=5)

    fig = plt.figure(dpi=300, figsize=(fig_width / 2.54, fig_height / 2.54))
    ax = fig.add_subplot()

    # IAM for PTC #######################################################
    ax.plot(angles, ones(shape=angles.shape[0]), label='PTC - $K_T$',
            color='tab:blue',
            marker='s', ms=ms,
            lw=lw)
    ax.plot(angles, ptc_kl_f(angles), label='PTC - $K_L$',
            color='tab:blue',
            marker='^', ms=ms,
            ls='dashed', lw=lw)

    # IAM for LFC ########################################################
    ax.plot(angles, lfc_kt_f(angles), label='LFC - $K_T$',
            color='tab:orange',
            marker='s', ms=ms,
            lw=lw)
    ax.plot(angles, lfc_kl_f(angles), label='LFC - $K_L$',
            color='tab:orange',
            marker='^', ms=ms,
            ls='dashed', lw=lw)
    ###################################################################################

    ax.set_xlabel(r'$\theta_T, \theta_L$', fontsize=font_size)
    ax.set_ylabel(r'IAM', fontsize=font_size)
    ax.xaxis.set_major_formatter(EngFormatter(unit=r"$^{\circ}$", sep=""))

    ax.set_xticks(arange(start=0, stop=95, step=10))

    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax.legend(fontsize=font_size)
    ax.grid(alpha=0.5)

    plt.tight_layout(pad=0)
    plt.savefig('ptc_vs_lfc_iam_comparison.pdf')
    plt.show()
