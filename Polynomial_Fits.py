

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ZScaleInterval
from scipy.optimize import curve_fit
import scipy.stats as stats


source = 4 #do not change to 1! 
degree = 5

frequencies_source = [908037109.375, 952341796.875, 996646484.375, 1043458984.375, 1092779296.875, 1144607421.875, 1317228515.625, 1381177734.375, 1448052734.375, 1519943359.375, 1593923828.125, 1656201171.875]
flux_densities_source1 = [0.0104139, 0.013042, 0.0116122, 0.00973644, 0.00934067, 0.00896728, 0.00821659, 0.00801487, 0.00744964, 0.00709993, 0.0064038, 0.00629778]
flux_densities_source2 = [0.00526457, 0.00561876, 0.00576872, 0.00344533, 0.00308578, 0.00276929, 0.00309652, 0.00331458, 0.00302509, 0.00267429, 0.00188835, 0.00195252 ]
flux_densities_source3 = [0.00338948, 0.004893, 0.00367111, 0.00313642, 0.00262961, 0.00201451, 0.00104168, 0.00119817, 0.000969179, 0.00082944, 0.000806129, 0.000864887 ]
flux_densities_source4 = [0.00686867, 0.00782422, 0.00739618, 0.00698157, 0.00587714, 0.00507117, 0.00440633, 0.00494582, 0.00503633, 0.0052514, 0.00528332, 0.00531395]
flux_densities_source5 = [0.0164708, 0.014907, 0.0140903, 0.011875, 0.0122675, 0.0111966, 0.00983453, 0.0103902, 0.0101012, 0.00989463, 0.00924476, 0.00904366]

Flux = [flux_densities_source1, flux_densities_source2, flux_densities_source3, flux_densities_source4, flux_densities_source5]

polynomial_coefficients = np.polyfit(frequencies_source, Flux[source-1], degree)
polynomial_fit = np.poly1d(polynomial_coefficients)

frequencies_fit = np.linspace(min(frequencies_source), max(frequencies_source), 500)
flux_fit = polynomial_fit(frequencies_fit)

a,b, c, d,e, f = polynomial_coefficients

expected_flux = polynomial_fit(frequencies_source)
observed_flux = np.array(Flux[source-1])
chi_squared = np.sum((observed_flux - expected_flux) ** 2)

n_data_points = len(frequencies_source)
num_parameters = degree +1  #polynomial degree + 1 for the constant term
degrees_of_freedom = n_data_points - num_parameters

p_value = 1 - stats.chi2.cdf(chi_squared, degrees_of_freedom)


print(f"Chi-squared: {chi_squared:.3e}")
print(f"Degrees of Freedom: {degrees_of_freedom}")
print(f"P-value: {p_value:.3e}")



plt.scatter(frequencies_source, Flux[source-1], label='Data')
plt.plot(frequencies_fit, flux_fit, 'r-', label=f'Fit: y = {a:.3e}x⁵ + {b:.3e}x⁴ + {c:.3e}x³ + {d:.3e}x² + {e:.3e}x + {f:.3e}')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Flux Density (Jy)')
plt.title('Flux Density vs. Frequency with Polynomial Fit of Point Source %i' %source)


plt.legend(loc='upper center', bbox_to_anchor=(0.4, -0.15), fancybox=True, shadow=True, ncol=1, fontsize='small', handlelength=2)


plt.tight_layout(rect=[0, 0.1, 1, 1])
plt.savefig('Flux_Density_vs_Freq_source%i.png' %source)
plt.show()
