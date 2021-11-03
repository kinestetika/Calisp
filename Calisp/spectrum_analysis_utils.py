import numpy as np
from scipy.optimize import minimize_scalar,OptimizeResult
from Calisp import element_count_and_mass_utils as utils

ELEMENT_ROW_INDEX = 0
ISOTOPE_COLUMN_INDEX = 1

MODEL_NEUTRON_ABUNDANCE = 0
MODEL_FFT_SS = 1
MODEL_FFT_NORM = 2
MODEL_CLUMPY_CARBON = 3
MODEL_WOBBLINESSS = 4

MASS_SHIFT_VARIATION_TOLERANCE = 1e-5
NEUTRON_MASS_SHIFT_TOLERANCE = 0.002

ISOTOPE_MATRIX = np.array([
        [0.988943414833479, 0.011056585166521, 0.0,         0.0, 0.0,    0.0, 0.0], # C
        [0.996323567,       0.003676433,       0.0,         0.0, 0.0,    0.0, 0.0], # N
        [0.997574195,       0.00038,           0.002045805, 0.0, 0.0,    0.0, 0.0], # O
        [0.99988,           0.00012,           0.0,         0.0, 0.0,    0.0, 0.0], # H
        [0.9493,            0.0076,            0.0429,      0.0, 0.0002, 0.0, 0.0]  # S
    ], dtype=np.float32)

NATURAL_ABUNDANCES = np.array([
        [0.988943414833479, 0.011056585166521, 0.0,         0.0, 0.0,    0.0, 0.0], # C
        [0.996323567,       0.003676433,       0.0,         0.0, 0.0,    0.0, 0.0], # N
        [0.997574195,       0.00038,           0.002045805, 0.0, 0.0,    0.0, 0.0], # O
        [0.99988,           0.00012,           0.0,         0.0, 0.0,    0.0, 0.0], # H
        [0.9493,            0.0076,            0.0429,      0.0, 0.0002, 0.0, 0.0]  # S
    ], dtype=np.float32)

# VPDB standard 13C/12C = 0.0111802 in Isodat software
# https://www.webelements.com/sulfur/isotopes.html
# see also http://iupac.org/publications/pac/pdf/2003/pdf/7506x0683.pdf


def compute_relative_neutron_abundance_of_non_target_isotopes(element_counts, matrix=NATURAL_ABUNDANCES):
    rna_untargeted = np.float32(0)

    for row in range(len(matrix)):
        for column in range(len(matrix[0])):
            if row == ELEMENT_ROW_INDEX and column == ISOTOPE_COLUMN_INDEX:
                continue
            rna_untargeted += column * element_counts[row] * matrix[row][column]
    return rna_untargeted


def compute_relative_neutron_abundance(matrix, peaks, element_counts):
    rna = np.sum(peaks * range(len(peaks))) / np.sum(peaks)
    rna -= compute_relative_neutron_abundance_of_non_target_isotopes(element_counts)
    rna /= element_counts[ELEMENT_ROW_INDEX]
    rna /= ISOTOPE_COLUMN_INDEX
    return rna


def standard_ratio():
    return ISOTOPE_MATRIX[ELEMENT_ROW_INDEX][ISOTOPE_COLUMN_INDEX] / ISOTOPE_MATRIX[ELEMENT_ROW_INDEX][0]


def __fft(element_counts:np.ndarray, matrix:np.ndarray, modeled_peaks:np.ndarray):
    fft_vector_size = len(modeled_peaks)
    transforms_1 = np.empty((len(matrix), fft_vector_size), dtype=np.complex64)
    for i in range(len(matrix)): # iterate over rows (elements)
        v = np.empty(len(matrix[0]), dtype = np.float32)
        for j in range(len(matrix[0])): # iterate over columns (isotopes)
            v[j] = matrix[i][j]
        V: np.ndarray = np.fft.fft(v, n=fft_vector_size) # norm{“backward”, “ortho”, “forward”}, optional
        for j in range(fft_vector_size): # iterate over spectrum
            V[j] = pow(V[j], element_counts[i])
        transforms_1[i] = V
    transforms_2 = np.empty(fft_vector_size, dtype=np.complex64)
    for j in range(fft_vector_size): # iterate over spectrum
        W = complex(1,0)
        for i in range(len(matrix)): # iterate over rows (elements):
            W *= transforms_1[i][j]
        transforms_2[j] = W
    X: np.ndarray = np.fft.ifft(transforms_2) # inverse FFT
    total_intensity = np.float32(0)
    for j in range(fft_vector_size): # iterate over spectrum
        total_intensity += X[j].real
    for j in range(fft_vector_size): # iterate over spectrum
        modeled_peaks[j] = X[j].real / total_intensity


def __fft_vector_len(peak_count, element_counts:np.ndarray):
    fft_vector_len = 128
    if peak_count <= 48:
        fft_vector_len = 64
    if peak_count <= 24:
        fft_vector_len = 32
    if peak_count <= 12:
        fft_vector_len = 16
    if peak_count <= 4:
        fft_vector_len = 8
    modeled_peaks = np.zeros(fft_vector_len, dtype=np.float32)
    __fft(element_counts, ISOTOPE_MATRIX, modeled_peaks)
    while modeled_peaks[-1] > 1e-6:
        fft_vector_len *= 2
        modeled_peaks = np.zeros(fft_vector_len, dtype=np.float32)
        __fft(element_counts, ISOTOPE_MATRIX, modeled_peaks)
        if fft_vector_len >= 256:
            break
    return fft_vector_len


def __fft_fitting_function(rna, element_counts, matrix, experimental_peaks, modeled_peaks):
    matrix[ELEMENT_ROW_INDEX][ISOTOPE_COLUMN_INDEX] = rna
    matrix[ELEMENT_ROW_INDEX][0] = 1 - np.sum(matrix[ELEMENT_ROW_INDEX][1:])
    __fft(element_counts, matrix, modeled_peaks)
    return distance_ss(experimental_peaks, modeled_peaks)


# def __fft_fitting_function_clumpy_carbon(current_ratio, i, element_counts, matrix, ratios,
#                                          fft_vector_size, experimental_peaks):
#     ratios[i] = current_ratio
#     matrix[ELEMENT_ROW_INDEX][:i] = ratios[:i] * matrix[ELEMENT_ROW_INDEX][0] / sum(ratios[:i])
#     __fft(element_counts, matrix, fft_vector_size)
#     return distance_ss(experimental_peaks, model_peaks, i + 1)


def distance_ss(peaks_a:np.ndarray, peaks_b:np.ndarray, max_peak_count=9999): ##assumes normalized spectra
    shared_peak_count = min(len(peaks_a), len(peaks_b), max_peak_count)
    if not shared_peak_count:
        return np.float32(0)

    peaks_a = peaks_a[:shared_peak_count]
    peaks_b = peaks_b[:shared_peak_count]

    total_a = np.sum(peaks_a)
    total_b = np.sum(peaks_b)
    if not total_a and not total_b:
        return 0
    elif not total_a:
        return 9999
    elif not total_b:
        return 9999
    peaks_a = peaks_a / total_a
    peaks_b = peaks_b / total_b
    sum_of_squares = np.float32(0)
    for i in range(shared_peak_count):
        sum_of_squares += (peaks_a[i] - peaks_b[i]) ** 2
    return sum_of_squares


def fit_relative_neutron_abundance(spectrum:{}, experimental_peaks:np.ndarray, element_counts:np.ndarray,
                                   matrix:np.ndarray = ISOTOPE_MATRIX):
    if not sum(element_counts):
        return {}
    rna = compute_relative_neutron_abundance(matrix, experimental_peaks, element_counts)
    ratio = rna / (1 - rna)
    spectrum['ratio_na'] = ratio
    return spectrum


def fit_fft(spectrum:{}, experimental_peaks:np.ndarray, element_counts:np.ndarray, matrix:np.ndarray = ISOTOPE_MATRIX):
    if not sum(element_counts):
        return {}
    fft_vector_size = __fft_vector_len(len(experimental_peaks), element_counts)
    modeled_peaks = np.zeros(fft_vector_size, dtype=np.float32)
    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalar.html
    fit:OptimizeResult = minimize_scalar(__fft_fitting_function, bounds=(0.0001, 0.15), method='Bounded',
                                         args=(element_counts, matrix, experimental_peaks, modeled_peaks), options={'maxiter':40})
    spectrum['ratio_fft'] = fit['x']
    spectrum['error_fft'] = fit['fun']
    return spectrum


def fit_clumpy_carbon(spectrum:{}, experimental_peaks:np.ndarray, matrix:np.ndarray, element_counts:np.ndarray):
    if not sum(element_counts):
        return {}
    fft_vector_size = __fft_vector_len(len(experimental_peaks), element_counts)
    matrix[ELEMENT_ROW_INDEX] = np.array([1,0,0,0,0,0,0], dtype=np.float32)
    ratios = np.array([1,0,0,0,0,0,0], dtype=np.float32)
    peak_count = min(6, len(experimental_peaks))
    modeled_peaks = np.zeros(fft_vector_size, dtype=np.float32)
    for i in range(1, peak_count):
        fit: OptimizeResult = minimize_scalar(__fft_fitting_function, bounds=(0.0001, 0.15), method='Bounded',
                                              args=(element_counts, matrix, fft_vector_size, experimental_peaks),
                                              options={'maxiter':40})
        #ratio=0
        #for i2 in range(1,i):
        #     ratio += i2 * matrix[ELEMENT_ROW_INDEX][i2]
        #ratio /= matrix[ELEMENT_ROW_INDEX][0]


def compute_spacing_and_irregularity(spectrum:{}, masses, charge):
    if not len(masses):
        return
    peak_mass_spacings = np.empty(len(masses)-1, dtype=np.float32)
    for i in range(1,len(masses)):
        peak_mass_spacings[i-1] = masses[i] - masses[i-1]
    spacing_count = len(peak_mass_spacings)
    if spacing_count % 2 == 0:
        median_peak_spacing = (peak_mass_spacings[int(spacing_count/2)] +
                               peak_mass_spacings[int(spacing_count/2)-1])/2
    else:
        median_peak_spacing = peak_mass_spacings[int(spacing_count/2)]
    mass_irregularity = np.float32(0)
    for i in range(1,len(masses)):
        spacing = masses[i] - masses[i-1]
        mass_irregularity += (spacing - median_peak_spacing)**2
    median_peak_spacing *= charge
    mass_irregularity *= charge
    mass_irregularity /= len(masses)
    is_wobbly = mass_irregularity > MASS_SHIFT_VARIATION_TOLERANCE and \
                abs(median_peak_spacing - utils.NEUTRON_MASS_SHIFT) < NEUTRON_MASS_SHIFT_TOLERANCE
    spectrum['spectrum_median_peak_spacing'] = median_peak_spacing
    spectrum['spectrum_mass_irregularity'] = mass_irregularity
    spectrum['flag_spectrum_is_wobbly'] = is_wobbly
    return spectrum
