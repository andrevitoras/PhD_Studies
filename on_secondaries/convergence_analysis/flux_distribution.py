import matplotlib.pyplot as plt
from matplotlib import ticker
from numpy import meshgrid, ones

from on_secondaries.convergence_analysis.base_models import *
from on_secondaries.convergence_analysis.classes_and_functions import *

simulation_files_path = Path(Path.cwd(), 'simulation_files', 'flux', 'test_cases')


trace_options = Trace(rays=1e6, seed=123, cpus=8,
                      simulate=True, optical_errors=True, sunshape=True)

absorber_data = simulate_optic(file_name='optic', file_path=simulation_files_path,
                               theta_t=0., theta_l=0., dni=1000.,
                               lfc=lfc_1, x_bins=150, y_bins=150,
                               source=sun_source, optics=optical_properties,
                               trace=trace_options,
                               soltrace_version=2012, force_sim=False)


x_values = absorber_data.xbin_values
y_values = absorber_data.ybin_values
r_values = absorber_data.rbin_values
X, Y = meshgrid(x_values, y_values)

fig1 = plt.figure(dpi=300)

ax1 = fig1.add_subplot()
c_map = ax1.contourf(X, Y, absorber_data.flux_map, cmap="turbo")

cbar = fig1.colorbar(c_map, ax=ax1)
plt.xlabel('$x$ [m]')
plt.ylabel('$y$ [m]')
plt.title('Flux intensity map')
cbar.set_label('$\mathrm{kW/m^2}$')

plt.show()

fig2 = plt.figure(dpi=300, figsize=(12, 5))

# a polar plot for the circumferential flux distribution
# ax2 = fig2.add_subplot(1, 2, 1, projection='polar')
# plt.polar(r_values * pi/180, absorber_data.flux_x_distribution)
# ax2.set_theta_zero_location('S')

# A linear plot for the circumferential flux distribution
ax2 = fig2.add_subplot(1, 2, 1)
plt.plot(r_values, absorber_data.flux_x_distribution)
ax2.hlines(y=absorber_data.flux_mean, xmin=r_values.min(), xmax=r_values.max(),
           color='red', label='Average flux intensity')
plt.xlabel(r'$\theta$ [degrees]')
plt.ylabel('Flux intensity [$\mathrm{kW/m^2}$]')
plt.legend()

ax2 = fig2.add_subplot(1, 2, 2)
plt.plot(y_values, absorber_data.flux_y_distribution)
ax2.hlines(y=absorber_data.flux_mean, xmin=y_values.min(), xmax=y_values.max(),
           color='red')

plt.xlabel(r'$y$ [m]')
plt.show()


