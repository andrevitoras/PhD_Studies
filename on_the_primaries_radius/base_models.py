from numpy import sqrt

from on_the_primaries_radius.classes_and_functions import *


#################################################################################
# Study cases definition ########################################################
# Geometries

lfc_1 = lfr_geometry(name='lfc_1',
                     mirror_width=750,
                     total_width=16560,
                     nbr_mirrors=16,
                     rec_width=340,
                     rec_center=array([0, 7900.]))

lfc_2 = lfr_geometry(name='lfc_2',
                     mirror_width=750,
                     total_width=16560,
                     nbr_mirrors=20,
                     rec_width=340,
                     rec_center=array([0, 7900.]))
################################################################################

# Effective sources ###############################################
# 4.65 mrad half-width pillbox without optical errors
ES1 = RadialSource(profile='p', size=4.65e-3, name='ES1')

# 2.8 mrad standard deviation Gaussian sunshape
# with 5.0 mrad Gaussian optical errors
g_size = sqrt(2.8**2 + 5.0**2)
ES2 = RadialSource(profile='g', size=g_size * 1e-3, name='ES2')
###################################################################


# Site weather data #############################################################
emsp_lat, emsp_long = 38.5, -8.0
evora = SiteData(name='Evora', latitude=emsp_lat, longitude=emsp_long)

aswan = SiteData(name='Aswan', latitude=24.1, longitude=32.9,
                 start_year=2005, end_year=2015)


if __name__ == '__main__':
    evora.tmy_data.to_csv('evora_tmy.csv')
    aswan.tmy_data.to_csv('aswan_tmy.csv')

    print(evora.dni_sum / 1000)
    print(aswan.dni_sum / 1000)
