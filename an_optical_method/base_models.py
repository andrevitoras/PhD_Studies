from solaren.scopy.sunlight import SiteData
from an_optical_method.classes_and_functions import *

lfc_1 = uniform_geometry(name='LFC_1',
                         number_of_mirrors=16,
                         mirror_width=750., mirror_shift=1054.,
                         receiver_height=7400.,
                         radius_design='zenithal',
                         length=30e3)


lfc_2 = uniform_geometry(name='LFC_2',
                         number_of_mirrors=11,
                         mirror_width=250., mirror_shift=275.,
                         receiver_height=4000.,
                         radius_design='flat',
                         length=30e3)

ES1 = EffectiveSource(name='ES1',
                      sun=SunSettings(sunshape=None, size=0.),
                      optical=OpticalSettings(slope_error=0., specular_error=0.))

ES2 = EffectiveSource(name='ES2',
                      sun=SunSettings(sunshape='p', size=4.65e-3),
                      optical=OpticalSettings(slope_error=0., specular_error=0.))

ES3 = EffectiveSource(name='ES3',
                      sun=SunSettings(sunshape='p', size=4.65e-3),
                      optical=OpticalSettings(slope_error=2e-3, specular_error=3.e-3))

ES4 = EffectiveSource(name='ES4',
                      sun=SunSettings(sunshape='g', size=2.8e-3),
                      optical=OpticalSettings(slope_error=2e-3, specular_error=3.e-3))


ES5 = EffectiveSource(name='ES5',
                      sun=SunSettings(sunshape='b', size=0.025),
                      optical=OpticalSettings(slope_error=2e-3, specular_error=3.e-3))


lodwar = SiteData(name='Lodwar (Kenya)', latitude=3.12, longitude=35.60)
denbel = SiteData(name='Denbel (Ethiopia)', latitude=9.85, longitude=42.78)
najran = SiteData(name='Najran (Saudi Arabia)', latitude=17.57, longitude=44.23)
aswan = SiteData(name='Aswan (Egypt)', latitude=24.1, longitude=32.9)
quarzazate = SiteData(name='Quarzazate (Morocco)', latitude=31.01, longitude=-6.86)
evora = SiteData(name='Evora (Portugal)', latitude=38.5, longitude=-8.0)

geometries = [lfc_1, lfc_2]
effective_sources = [ES1, ES2, ES3, ES4, ES5]
selected_sites = [lodwar, denbel, najran, aswan, quarzazate, evora]

if __name__ == '__main__':

    [print(f'The RMS width of {es.name} is {round(es.rms_width*1000, 2)} mrad!')
     for es in effective_sources]

    [print(f'The DNI sum for {loc.name} is {round(loc.dni_sum/1000., 2)} kW/m2.year!')
     for loc in selected_sites]
