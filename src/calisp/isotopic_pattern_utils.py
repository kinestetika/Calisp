from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar, OptimizeResult

from calisp import element_count_and_mass_utils

ELEMENT_ROW_INDEX = 0
ISOTOPE_COLUMN_INDEX = 1

MASS_SHIFT_VARIATION_TOLERANCE = 1e-5
NEUTRON_MASS_SHIFT_TOLERANCE = 0.002

CLUMPY_CARBON_BOUNDS = [0.1, 0.1, 0.05, 0.0375, 0.02, 0.02, 0.02]

DEFAULT_MATRIX_FILE = Path(__file__).parent / 'isotope_matrix.txt'
ISOTOPE_MATRIX = None #isotopic_pattern_utils.load_isotope_matrix(DEFAULT_MATRIX_FILE)
NATURAL_ABUNDANCES = None #isotopic_pattern_utils.load_isotope_matrix(DEFAULT_MATRIX_FILE)


def load_isotope_matrix(matrix_file):
    if not matrix_file:
        matrix_file = DEFAULT_MATRIX_FILE
    matrix = np.zeros((5, 7), dtype=np.float32)
    if not matrix_file.exists():
        raise Exception(f'Isotope matrix file {matrix_file} not found!')
    with open(matrix_file) as matrix_reader:
        row = 0
        for line in matrix_reader:
            if line.startswith('#'):
                continue
            if (comment_index := line.find('#')) >= 0:
                line = line[0:comment_index]
            line_data = line.split()
            for column in range(0, min(len(line_data), 6)):
                matrix[row, column] = float(line_data[column])
            row += 1
    # test if matrix is sane:
    for row in range(5):
        for column in range(2):
            if not matrix[row, column]:
                raise Exception (f"Isotope matrix has zero at row {row}, column {column}.")
    return matrix


def print_matrix(isotope_column_index=-1, element_row_index=-1, matrix=None):
    if not matrix:
        matrix = ISOTOPE_MATRIX
    elements = 'CNOHS'
    print('   isotope matrix')
    for row_id in range(len(matrix)):
        print('   ', end='')
        end = '* ' if row_id == element_row_index else '  '
        print(elements[row_id], end=end)
        row = matrix[row_id]
        for column_id in range(len(row)):
            end = '* ' if row_id == element_row_index and column_id == isotope_column_index else '  '
            print(f'{row[column_id]:.5f}', end=end)
        print()


def compute_relative_neutron_abundance_of_non_target_isotopes(element_counts, matrix: np.ndarray):
    rna_untargeted = np.float32(0)
    for row in range(len(matrix)):
        for column in range(len(matrix[0])):
            if row == ELEMENT_ROW_INDEX and column == ISOTOPE_COLUMN_INDEX:
                continue
            rna_untargeted += column * element_counts[row] * matrix[row][column]
    return rna_untargeted


def compute_relative_neutron_abundance(peaks, element_counts, matrix: np.ndarray):
    rna = np.sum(peaks * range(len(peaks))) / np.sum(peaks)
    rna -= compute_relative_neutron_abundance_of_non_target_isotopes(element_counts, matrix)
    rna /= element_counts[ELEMENT_ROW_INDEX]
    rna /= ISOTOPE_COLUMN_INDEX
    return rna


def standard_ratio(matrix):
    return matrix[ELEMENT_ROW_INDEX][ISOTOPE_COLUMN_INDEX] / matrix[ELEMENT_ROW_INDEX][0]


def __fft(element_counts: np.ndarray, matrix: np.ndarray, modeled_peaks: np.ndarray):
    fft_vector_size = len(modeled_peaks)
    transforms_1 = np.empty((len(matrix), fft_vector_size), dtype=np.complex64)
    for i in range(len(matrix)):  # iterate over rows (elements)
        v = np.empty(len(matrix[0]), dtype=np.float32)
        for j in range(len(matrix[0])):  # iterate over columns (isotopes)
            v[j] = matrix[i][j]
        v_v: np.ndarray = np.fft.fft(v, n=fft_vector_size)  # norm{“backward”, “ortho”, “forward”}, optional
        for j in range(fft_vector_size):  # iterate over spectrum
            v_v[j] = pow(v_v[j], element_counts[i])
        transforms_1[i] = v_v
    transforms_2 = np.empty(fft_vector_size, dtype=np.complex64)
    for j in range(fft_vector_size):  # iterate over spectrum
        v_w = complex(1, 0)
        for i in range(len(matrix)):  # iterate over rows (elements):
            v_w *= transforms_1[i][j]
        transforms_2[j] = v_w
    v_x: np.ndarray = np.fft.ifft(transforms_2)  # inverse FFT
    total_intensity = np.float32(0)
    for j in range(fft_vector_size):  # iterate over spectrum
        total_intensity += v_x[j].real
    for j in range(fft_vector_size):  # iterate over spectrum
        modeled_peaks[j] = v_x[j].real / total_intensity


def __fft_vector_len(peak_count, element_counts: np.ndarray, matrix: np.ndarray):
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
    __fft(element_counts, matrix, modeled_peaks)
    while modeled_peaks[-1] > 1e-6:
        fft_vector_len *= 2
        modeled_peaks = np.zeros(fft_vector_len, dtype=np.float32)
        __fft(element_counts, matrix, modeled_peaks)
        if fft_vector_len >= 256:
            break
    return fft_vector_len


def __fft_fitting_function(rna, element_counts, matrix, experimental_peaks, modeled_peaks):
    matrix[ELEMENT_ROW_INDEX][ISOTOPE_COLUMN_INDEX] = rna
    matrix[ELEMENT_ROW_INDEX][0] = 1 - np.sum(matrix[ELEMENT_ROW_INDEX][1:])
    __fft(element_counts, matrix, modeled_peaks)
    return distance_ss(experimental_peaks, modeled_peaks)


def __fft_fitting_function_clumpy_carbon(rna, element_counts, matrix, experimental_peaks, modeled_peaks, i):
    matrix[ELEMENT_ROW_INDEX][i] = rna
    matrix[ELEMENT_ROW_INDEX][:i] = matrix[ELEMENT_ROW_INDEX][:i] * (1-rna) / sum(matrix[ELEMENT_ROW_INDEX][:i])
    __fft(element_counts, matrix, modeled_peaks)
    return distance_ss(experimental_peaks, modeled_peaks, i+1)


def distance_ss(peaks_a: np.ndarray, peaks_b: np.ndarray, max_peak_count=9999):  # assumes normalized spectra
    shared_peak_count = min(len(peaks_a), len(peaks_b), max_peak_count)
    if not shared_peak_count:
        return 9999

    peaks_a = peaks_a[:shared_peak_count]
    peaks_b = peaks_b[:shared_peak_count]

    total_a = np.sum(peaks_a)
    total_b = np.sum(peaks_b)
    if not total_a and not total_b:
        return 9999
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


def fit_relative_neutron_abundance(spectrum: {}, experimental_peaks: np.ndarray, element_counts: np.ndarray,
                                   matrix: np.ndarray):
    if not sum(element_counts):
        return {}
    rna = compute_relative_neutron_abundance(experimental_peaks, element_counts, matrix)
    ratio = rna / (1 - rna)
    spectrum['ratio_na'] = ratio
    return spectrum


def fit_fft(pattern: {}, experimental_peaks: np.ndarray, element_counts: np.ndarray,
            matrix: np.ndarray):
    if not sum(element_counts):
        return {}
    fft_vector_size = __fft_vector_len(len(experimental_peaks), element_counts, matrix)
    modeled_peaks = np.zeros(fft_vector_size, dtype=np.float32)
    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize_scalar.html
    fit: OptimizeResult = minimize_scalar(__fft_fitting_function, bounds=(0, 0.15), method='Bounded',
                                          args=(element_counts, matrix, experimental_peaks, modeled_peaks),
                                          options={'maxiter': 40})
    pattern['ratio_fft'] = fit['x']
    pattern['error_fft'] = fit['fun']
    return pattern


def fit_clumpy_carbon(pattern: {}, experimental_peaks: np.ndarray, element_counts: np.ndarray,
                      matrix: np.ndarray):
    if not sum(element_counts):
        return {}
    fft_vector_size = __fft_vector_len(len(experimental_peaks), element_counts, matrix)
    matrix[ELEMENT_ROW_INDEX] = np.array([1, 0, 0, 0, 0, 0, 0], dtype=np.float32)
    index_len = min(6, len(experimental_peaks)-1)
    modeled_peaks = np.zeros(fft_vector_size, dtype=np.float32)
    for i in range(1, index_len):
        fit: OptimizeResult = minimize_scalar(__fft_fitting_function_clumpy_carbon, bounds=(0, CLUMPY_CARBON_BOUNDS[i]),
                                              method='Bounded', args=(element_counts, matrix, experimental_peaks,
                                                                      modeled_peaks, i),
                                              options={'maxiter': 40})
        pattern['error_clumpy'] = fit['fun']
    ratios = matrix[ELEMENT_ROW_INDEX][1:]*range(1, len(matrix[ELEMENT_ROW_INDEX]))
    ratios /= sum(ratios)
    for i in range(index_len):
        pattern[f'c{i + 1}'] = ratios[i]
    return pattern


def compute_spacing_and_irregularity(pattern: {}, masses, charge):
    if not len(masses):
        return
    peak_mass_spacings = np.empty(len(masses)-1, dtype=np.float32)
    for i in range(1, len(masses)):
        peak_mass_spacings[i-1] = masses[i] - masses[i-1]
    spacing_count = len(peak_mass_spacings)
    if spacing_count % 2 == 0:
        median_peak_spacing = (peak_mass_spacings[int(spacing_count/2)] +
                               peak_mass_spacings[int(spacing_count/2)-1])/2
    else:
        median_peak_spacing = peak_mass_spacings[int(spacing_count/2)]
    mass_irregularity = np.float32(0)
    for i in range(1, len(masses)):
        spacing = masses[i] - masses[i-1]
        mass_irregularity += (spacing - median_peak_spacing)**2
    median_peak_spacing *= charge
    mass_irregularity *= charge
    mass_irregularity /= len(masses)
    is_wobbly = mass_irregularity > MASS_SHIFT_VARIATION_TOLERANCE and \
        abs(median_peak_spacing - element_count_and_mass_utils.NEUTRON_MASS_SHIFT) < NEUTRON_MASS_SHIFT_TOLERANCE
    pattern['pattern_median_peak_spacing'] = median_peak_spacing
    pattern['pattern_mass_irregularity'] = mass_irregularity
    pattern['flag_pattern_is_wobbly'] = is_wobbly
    return pattern
