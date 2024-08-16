# Written by Lynn-Joy and Abigail

# importing necessary libraries
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ZScaleInterval

with fits.open('/home/abigail/Hons/OT2/Assignment_2/Abell_2744_aFix_pol_I_15arcsec_fcube_cor.fits') as hdul:
    data_cube = hdul[0].data
    header = hdul[0].header

print(np.shape(data_cube))
print(header)

for i in range(1, 13):
    
    slice_data = data_cube[0][i-1][:][:]
    slice_header = header.copy()
    
    for j in range(1, 13):
        if j != i:
            slice_header[f'FREQ{j:04d}'] = 'Not Available'
    
    frequency = header.get(f'FREQ{i:04d}', 'Keyword not found')
    slice_header[f'FREQ{i:04d}'] = frequency

    # saving the slice data to a new FITS file
    hdu = fits.PrimaryHDU(data=slice_data, header=slice_header)
    hdu.writeto(f'slice_{i}.fits', overwrite=True)
    
    freq = f"{frequency:.3e}"
    print()
    print(f'FREQ{i:04d} (frequency of slice {i}): {frequency}')
    
    plt.imshow(slice_data, vmin=ZScaleInterval().get_limits(slice_data)[0], vmax=ZScaleInterval().get_limits(slice_data)[1], origin='lower', cmap='magma')
    plt.colorbar()
    plt.title(f'Slice {i} at frequency {freq} Hz')
    
    # saving the image before showing it (hashtaged so that it only saves if I need it to)
    #plt.savefig(f'slice_image_{i}.png')
    
    plt.show()
    
    plt.close()
    

hdul.close()

freq_source = [
    908037109.375, 952341796.875, 996646484.375, 1043458984.375,
    1092779296.875, 1144607421.875, 1317228515.625, 1381177734.375,
    1448052734.375, 1519943359.375, 1593923828.125, 1656201171.875
]

flux_dens_source1 = [0.0104139, 0.013042, 0.0116122, 0.00973644, 0.00934067, 0.00896728, 0.00821659, 0.00801487, 0.00744964, 0.00709993, 0.0064038, 0.00629778]
flux_dens_source2 = [0.00526457, 0.00561876, 0.00576872, 0.00344533, 0.00308578, 0.00276929, 0.00309652, 0.00331458, 0.00302509, 0.00267429, 0.00188835, 0.00195252]
flux_dens_source3 = [0.00338948, 0.004893, 0.00367111, 0.00313642, 0.00262961, 0.00201451, 0.00104168, 0.00119817, 0.000969179, 0.00082944, 0.000806129, 0.000864887]
flux_dens_source4 = [0.00686867, 0.00782422, 0.00739618, 0.00698157, 0.00587714, 0.00507117, 0.00440633, 0.00494582, 0.00503633, 0.0052514, 0.00528332, 0.00531395]
flux_dens_source5 = [0.0164708, 0.014907, 0.0140903, 0.011875, 0.0122675, 0.0111966, 0.00983453, 0.0103902, 0.0101012, 0.00989463, 0.00924476, 0.00904366]

flux_dens = [flux_dens_source1, flux_dens_source2, flux_dens_source3, flux_dens_source4, flux_dens_source5]

for i in range(5):
    flux = flux_dens[i]
    flux1 = flux[0]
    flux2 = flux[-1]
    
    freq1 = freq_source[0]
    freq2 = freq_source[-1]
    
    spectral_index = np.log(flux2 / flux1) / np.log(freq2 / freq1)
    print(f"Spectral index of source {i+1}: {spectral_index:.4f}")