from pathlib import Path

import matplotlib.pyplot as plt
from numpy import array
from pandas import read_csv, DataFrame
from scipy.optimize import root

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\\usepackage{amsmath}\n \\usepackage{amssymb}'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'NewComputerModern10'

########################################################################################################################
# Global variables for this program ####################################################################################

figures_path = Path(Path.cwd(), 'figures')
figures_path.mkdir(parents=True, exist_ok=True)

reference_colors = {'PTC': 'tab:blue',
                    'CRS': 'tab:orange',
                    'LFC': 'tab:green'}


########################################################################################################################
########################################################################################################################

########################################################################################################################
# Functions for tis program ############################################################################################


def classify_technology(tech):
    if 'Trough' in tech:
        return 'PTC'
    elif 'Tower' in tech:
        return 'CRS'
    elif 'Fresnel' in tech:
        return 'LFC'
    else:
        return 'Other'


def read_csp_database(file_full_path):
    df = read_csv(file_full_path)
    df['Technology_Class'] = df['Technology'].apply(lambda x: classify_technology(str(x)))

    return df


def plot_ste_plants(
        database: DataFrame,
        figwidth=12, figheight=6.,
        font_size=9,
        fig_name='csp_technologies', fig_format='pdf'):

    df = database
    df['Year_Temporal'] = df['Year_construction_start'].fillna(df['Year_operational'])

    # Drop rows with missing years
    df_filtered = df.dropna(subset=['Year_Temporal'])

    # Filter out 'Other' category
    df_filtered = df_filtered[df_filtered['Technology_Class'] != 'Other']

    # Group by year and sum the capacity
    cumulative_capacity = df_filtered.groupby(
        ['Year_Temporal', 'Technology_Class'])['Capacity_MW'].sum().unstack().fillna(0)
    # Convert the year column to integer for better plotting
    cumulative_capacity.index = cumulative_capacity.index.astype(int)

    # Compute cumulative sum over time
    cumulative_capacity = cumulative_capacity.cumsum()

    tech_capacities = [cumulative_capacity[tech].values[-1] / 1000.
                       for tech in ['PTC', 'CRS', 'LFC']]

    print(f'Cumulative capacity is {array(tech_capacities).sum().round(1)} GW!')

    fig = plt.figure(dpi=300,
                     figsize=(figwidth / 2.54, figheight / 2.54))

    # the first plot ---------------------------------------------------------------------------------------------------
    ax = plt.subplot2grid(fig=fig, shape=(1, 3),
                          loc=(0, 0), colspan=2)

    [ax.plot(
        cumulative_capacity.index, cumulative_capacity[tech] * 1.e-3,
        color=reference_colors[tech],
        marker='o', ms=1.25, lw=0.75,
        label=tech)
        for tech in ['PTC', 'CRS', 'LFC']]

    ax.tick_params(axis='both', which='major', labelsize=font_size - 1)

    # ax.set_xlabel("Year", fontsize=font_size)
    ax.set_ylabel("Cumulative Capacity (GW)", fontsize=font_size)
    ax.set_title('(a)', fontsize=font_size)

    ax.legend(fontsize=font_size - 1)
    ax.grid(alpha=0.35)
    # ------------------------------------------------------------------------------------------------------------------

    # the second plot ------------------------------------------------------- ------------------------------------------
    ax = plt.subplot2grid(fig=fig, shape=(1, 3),
                          loc=(0, 2), colspan=1)

    operational_capacity = df_filtered[df_filtered['Status'] == 'Operational']

    tech_capacities = [operational_capacity[
                           operational_capacity['Technology_Class'] == tech]['Capacity_MW'].sum() / 1000.
                       for tech in ['PTC', 'CRS', 'LFC']]

    print(f'Operational capacity is {array(tech_capacities).sum().round(1)} GW!')

    ax.bar(x=[1, 2, 3],
           height=tech_capacities,
           color='#792459ff', alpha=1)

    ax.set_xticks([i + 1 for i, es in enumerate(tech_capacities)])
    ax.set_xticklabels([tech for i, tech in enumerate(['PTC', 'CRS', 'LFC'])])
    ax.grid(alpha=0.15)

    ax.tick_params(axis='both', which='major', labelsize=font_size - 1)

    ax.set_ylabel("Operational Capacity (GW)", fontsize=font_size)
    ax.set_title('(b)', fontsize=font_size)
    # ------------------------------------------------------------------------------------------------------------------

    plt.tight_layout(pad=0.5)
    plt.savefig(Path(figures_path, f'{fig_name}.{fig_format}'))
    plt.show()

    return cumulative_capacity


def plot_IEA_projections(figwidth=12, figheight=6.,
                         font_size=9,
                         fig_name='IEA_STE_capacity_vs_scenarios', fig_format='pdf'):

    """
    Data from IEA [1, pp. 299-311], tables with 'World electricity sector' data.

    [1] IEA, World Energy Outlook 2024, International Energy Agency (IEA), Paris, Tech. Rep., 2024.
        Available: https://www.iea.org/reports/world-energy-outlook-2024.

    :param figwidth:
    :param figheight:
    :param font_size:
    :param fig_name:
    :param fig_format:
    :return:
    """

    data_dic = {'years': [2023, 2030, 2035, 2040, 2050],
                'STEP': [7., 10., 21., 35., 68.],
                'APS': [7., 16., 57., 120., 230.],
                'NZE': [7., 35., 112., 226., 390.]
                }

    keys = ['STEP', 'APS', 'NZE']
    markers = ['o', 's', "^"]
    colors = ['#2c94abff', '#e16b09ff', '#25aa56ff']  # colors taken from the corresponding tables.

    fig = plt.figure(dpi=300,
                     figsize=(figwidth / 2.54, figheight / 2.54))
    ax = fig.add_subplot()

    [ax.plot(data_dic['years'], data_dic[k],
             label=k, lw=0.75, marker=m, ms=3, color=c)
     for k, m, c in zip(keys, markers, colors)]

    ax.legend(fontsize=font_size-1)
    ax.grid(alpha=0.25)

    ax.tick_params(axis='both', which='major', labelsize=font_size - 1)
    ax.set_ylabel('GW', fontsize=font_size)
    ax.set_xlabel('Year', fontsize=font_size)
    ax.set_title('Operating STE capacity', fontsize=font_size)

    plt.tight_layout(pad=0.5)
    plt.savefig(Path(figures_path, f'{fig_name}.{fig_format}'))
    plt.show()

    return None


def calculate_grown_rate(
        c0: float, cf: float,
        n: int, r0=0.2):

    r = root(lambda x: (cf/c0) - (1 + x)**n, x0=r0).x

    return r[0]


if __name__ == '__main__':
    database_file_name = 'csp-guru.csv'
    database_full_path = Path(Path.cwd(), database_file_name)
    csp_database = read_csp_database(file_full_path=database_full_path)

    plot_ste_plants(
        database=csp_database,
        figwidth=12, figheight=6.,
        font_size=9, fig_format='pdf')

    plot_IEA_projections(figwidth=8, figheight=7, fig_format='pdf')

    a = calculate_grown_rate(c0=7.0, cf=68, n=27)
    print(a * 100)
