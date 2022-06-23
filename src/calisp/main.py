import argparse
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pymzml
import time
from tqdm import tqdm

from calisp import isotopic_pattern_utils
from calisp import data_store
from calisp import element_count_and_mass_utils
from calisp.peptide_spectrum_match_files import PeptideSpectrumMatchFileReader

VERSION = '3.0.9'
START_TIME = time.monotonic()
LOG_TOPICS = set()
MASS_ACCURACY = 1e-5
COMPUTE_CLUMPS = False


def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def log(log_message, topic=''):
    if not topic or topic in LOG_TOPICS:
        print(f'{format_runtime()} {log_message}')


def parse_arguments():
    parser = argparse.ArgumentParser(description='src.py. (C) Marc Strous, Xiaoli Dong and Manuel Kleiner, 2018, 2021')
    parser.add_argument('--spectrum_file', required=True,  # type=argparse.FileType('r'),
                        help='[.mzML] file or folder with [.mzML] file(s).')
    parser.add_argument('--peptide_file', required=True,  # type=argparse.FileType('r'),
                        help='Search-engine-generated [.mzID] or [.target-peptide-spectrum-match] file or folder with'
                             ' these files.')
    parser.add_argument('--output_file', default='calisp-output',
                        help='The name of the folder the output gets written to. Default: calisp-output')
    parser.add_argument('--mass_accuracy', default=10, type=float,
                        help='The maximum mass difference between theoretical mass and experimental mass of a peptide')
    parser.add_argument('--bin_delimiter', default='_',
                        help='For metagenomic data, the delimiter that separates the bin ID from the protein ID '
                             '[default "_"]. Use "-" to ignore bins ID entirely.')
    parser.add_argument('--threads', default=4, type=int,
                        help='The number of (virtual) processors that calisp will use')
    parser.add_argument('--isotope', default='13C', choices=['13C', '14C', '15N', '17O', '18O', '2H', '3H', '33S',
                                                             '34S', '36S'],
                        help='The target isotope. Default: 13C')
    parser.add_argument('--compute_clumps', default=False, action='store_true',
                        help='To compute clumpiness of carbon assimilation. Only use when samples are labeled to'
                             'saturation. Estimation of clumpiness takes much additional time.')

    args = parser.parse_args()

    log(f'isotope:        {args.isotope} '
        f'[matrix{isotopic_pattern_utils.ELEMENT_ROW_INDEX, isotopic_pattern_utils.ISOTOPE_COLUMN_INDEX}]')
    log(f'peptide_file:   {args.peptide_file}')
    log(f'spectrum_file:  {args.spectrum_file}')
    log(f'output_file:    {args.output_file}')
    log(f'mass_accuracy:  {args.mass_accuracy} ppm')
    log(f'bin_delimiter:  {args.bin_delimiter}')
    log(f'threads:        {args.threads}')
    log(f'compute_clumps: {args.compute_clumps}')
    args.mass_accuracy /= 1e6
    return args


def prep_output_folder(output_folder):
    if os.path.exists(output_folder):
        log('Warning: may overwrite existing output files...')
        if os.path.isfile(output_folder):
            print(f'Expected folder at {output_folder}, found regular file, terminating now.')
            exit(1)
    else:
        os.mkdir(output_folder)


def flag_ambiguous_psms(my_psm_data: pd.DataFrame):
    for i in range(1, len(my_psm_data.index)):
        i1 = my_psm_data.index[i]
        i2 = my_psm_data.index[i-1]
        if my_psm_data.at[i1, 'psm_id'] == my_psm_data.at[i2, 'psm_id']:
            if my_psm_data.at[i1, 'psm_rank'] == my_psm_data.at[i2, 'psm_rank']:
                my_psm_data.at[i1, 'flag_psm_is_ambiguous'] = True
                my_psm_data.at[i2, 'flag_psm_is_ambiguous'] = True
            elif my_psm_data.at[i1, 'psm_rank'] > min(1, my_psm_data.at[i2, 'psm_rank']):
                my_psm_data.at[i1, 'flag_psm_is_ambiguous'] = True
            elif my_psm_data.at[i2, 'psm_rank'] > min(1, my_psm_data.at[i, 'psm_rank']):
                my_psm_data.at[i2, 'flag_psm_is_ambiguous'] = True


def set_precursor(precursors_list, my_psm_data, index):
    for precursor_dict in precursors_list:
        try:
            my_psm_data.at[index, 'psm_precursor_id'] = int(precursor_dict['precursor id'])
            my_psm_data.at[index, 'psm_precursor_mz'] = int(precursor_dict['mz'])
            return 1
        except KeyError:
            log('Warning: Unexpected definition of precursor dictionary:')
            for k in precursor_dict.keys():
                log(k, precursor_dict[k])
    return 0


def subsample_ms1_spectra(task):
    subsampled_patterns = []
    searches_done = set()
    for psm_index in task['psms'].index:
        # (1) subsample peptide isotopic patterns from precursor ms1 spectra
        precursor_id = task['psms'].at[psm_index, 'psm_precursor_id']
        if not precursor_id or task['psms'].at[psm_index, 'flag_peptide_mass_and_elements_undefined']:
            continue
        peptide_mass = task['psms'].at[psm_index, 'peptide_mass']
        precursor_i = task['ms1_spectra_hash'][precursor_id]
        initial_success = subsample_ms1_spectrum(precursor_i, psm_index, task, [], subsampled_patterns, peptide_mass,
                                                 searches_done)
        prev_success = initial_success
        # (success is a list of booleans, one boolean for each charge 1..5)
        for i in range(precursor_i+1, min(precursor_i+40, len(task['ms1_spectra_list']))):
            prev_success = subsample_ms1_spectrum(i, psm_index, task, prev_success, subsampled_patterns, peptide_mass,
                                                  searches_done)
            if not sum(prev_success):
                break
        prev_success = initial_success
        for i in range(precursor_i-1, max(0, precursor_i-40), -1):
            prev_success = subsample_ms1_spectrum(i, psm_index, task, prev_success, subsampled_patterns, peptide_mass,
                                                  searches_done)
            if not sum(prev_success):
                break
    for p in subsampled_patterns:
        # (2) assign isotopic patterns to the nearest psm:
        nearest_psm_index = p['psm_index']  # this is the psm for which the isotopic pattern was discovered
        shortest_distance = 999999
        if len(task['psms'].index) > 1:
            i_spe = task['ms1_spectra_hash'][p['pattern_precursor_id']]
            for psm_index in task['psms'].index:  # these are all the psms associated with the peptide
                i_psm = task['ms1_spectra_hash'][task['psms'].at[psm_index, 'psm_precursor_id']]
                distance_psm_to_pattern = (abs(i_psm - i_spe)+1) * (abs(task['psms'].at[psm_index, 'psm_charge']
                                                                         - p['pattern_charge']) + 1)
                if distance_psm_to_pattern < shortest_distance:
                    shortest_distance = distance_psm_to_pattern
                    nearest_psm_index = psm_index
            if nearest_psm_index != p['psm_index']:
                p['reassigned'] = True
        del p['psm_index']
        for key, value in task['psms'].loc[nearest_psm_index].to_dict().items():
            if key in p.keys():
                continue
            p[key] = value
        # (3) analyze isotopes:
        element_counts = np.array([p['C'], p['N'], p['O'], p['H'], p['S']], dtype=np.int16)
        intensities = np.empty(p['pattern_peak_count'], dtype=np.float32)
        masses = np.empty(p['pattern_peak_count'], dtype=np.float32)
        for i in range(p['pattern_peak_count']):
            intensities[i] = p[f'i{i}']
            masses[i] = p[f'm{i}']
        isotopic_pattern_utils.fit_relative_neutron_abundance(p, intensities, element_counts)
        isotopic_pattern_utils.fit_fft(p, intensities, element_counts)
        if COMPUTE_CLUMPS:
            isotopic_pattern_utils.fit_clumpy_carbon(p, intensities, element_counts)
        isotopic_pattern_utils.compute_spacing_and_irregularity(p, masses, p['pattern_charge'])
    return subsampled_patterns


def subsample_ms1_spectrum(i, psm_index, task, prev_success, subsampled_patterns, peptide_mass, searches_done):
    success = [False, False, False, False, False, False]
    if i in searches_done:
        return success
    searches_done.add(i)
    ms1_peaks = task['ms1_spectra_list'][i]
    for charge in range(1, 6):
        if len(prev_success) and not prev_success[charge]:
            continue
        expected_mz = peptide_mass / charge + element_count_and_mass_utils.PROTON_MASS
        # (1) find monoisotopic peak
        mono_isotopic_peak = []
        mip_p = 0
        for p in range(len(ms1_peaks)):
            peak = ms1_peaks[p]
            if abs(peak[0] - expected_mz) / expected_mz < MASS_ACCURACY:
                mono_isotopic_peak = peak
                mip_p = p
                break
            elif peak[0] > expected_mz:
                break
        if not len(mono_isotopic_peak):
            continue
        peak_moverz = [mono_isotopic_peak[0]]
        peak_intensities = [mono_isotopic_peak[1]]
        mass_shift = element_count_and_mass_utils.NEUTRON_MASS_SHIFT / charge

        # (2) assess presence of a peak at minus one position
        expected_mz = peak_moverz[0] - mass_shift
        peak_at_minus_one = False
        for q in range(mip_p-1, 0, -1):
            peak = ms1_peaks[q]
            if abs(peak[0] - expected_mz) / expected_mz < MASS_ACCURACY:
                peak_at_minus_one = True
            elif peak[0] < expected_mz:
                break

        # (3) find remainder of isotopic patterns
        for q in range(mip_p+1, len(ms1_peaks)):
            peak = ms1_peaks[q]
            # calculating distance from monoisotopic is slightly better than peak-to-peak
            expected_mz = peak_moverz[0] + mass_shift * len(peak_moverz)
            if abs(peak[0] - expected_mz) / expected_mz < MASS_ACCURACY:
                peak_moverz.append(peak[0])
                peak_intensities.append(peak[1])
            elif peak[0] > expected_mz:
                break

        if len(peak_moverz) >= 4:
            # normalize isotopic pattern
            peak_intensities = peak_intensities[:20]
            peak_moverz = peak_moverz[:20]
            total_intensity = sum(peak_intensities)
            peak_intensities /= total_intensity
            # assign isotopic pattern to the closest psm
            # bundle subsampled isotopic pattern as a dictionary
            subsampled_pattern = {'psm_index': psm_index,
                                   'pattern_charge': charge,
                                   'pattern_total_intensity': total_intensity,
                                   'pattern_precursor_id': task['ms1_id_list'][i],
                                   'pattern_peak_count': len(peak_intensities),
                                   'flag_pattern_is_contaminated': False,
                                   'flag_pattern_is_wobbly': False,
                                   'flag_peak_at_minus_one_pos': peak_at_minus_one,
                                   'reassigned': False}
            for p in range(len(peak_intensities)):
                subsampled_pattern[f'i{p}'] = peak_intensities[p]
                subsampled_pattern[f'm{p}'] = peak_moverz[p]
            subsampled_patterns.append(subsampled_pattern)
            success[charge] = True
    return success


def submit_for_overlap_detection(task):
    shared_patterns = task['patterns from same spectrum']
    for i in range(len(shared_patterns.index)):
        ii = shared_patterns.index[i]
        mii = set(shared_patterns.loc[ii,
                  'm0':f'm{shared_patterns.at[ii, "pattern_peak_count"] - 1}'].values)
        for j in range(i + 1, len(shared_patterns.index)):
            jj = shared_patterns.index[j]
            mjj = set(shared_patterns.loc[jj,
                      'm0':f'm{shared_patterns.at[jj, "pattern_peak_count"] - 1}'].values)
            if mjj.intersection(mii):
                return [ii, jj]
    return []


def main():
    log(f'This is calisp.py {VERSION}')
    args = parse_arguments()
    (isotopic_pattern_utils.ELEMENT_ROW_INDEX, isotopic_pattern_utils.ISOTOPE_COLUMN_INDEX) = \
        {'13C': (0, 1), '14C': (0, 2),
         '15N': (1, 1),
         '17O': (2, 1), '18O': (2, 2),
         '2H': (3, 1), '3H': (3, 2),
         '33S': (4, 1), '34S': (4, 2), '36S': (4, 4)}[args.isotope]
    global MASS_ACCURACY, COMPUTE_CLUMPS
    MASS_ACCURACY = args.mass_accuracy
    COMPUTE_CLUMPS = args.compute_clumps
    prep_output_folder(args.output_file)

    my_data_store = data_store.DataStore()
    my_data_store.set_scope(args.peptide_file, args.spectrum_file)

    log(f'Found {len(my_data_store.experiments)} peptide files and {len(my_data_store.ms_runs)} MS runs')
    for experiment_index in range(len(my_data_store.experiments)):

        # (1) read the "experiment" file with the psm and peptide definitions
        experiment_filename = my_data_store.experiments[experiment_index]
        log(f'Parsing ({experiment_index + 1}/{len(my_data_store.experiments)}) peptide file(s), '
            f'"{os.path.basename(experiment_filename)}"...')
        with PeptideSpectrumMatchFileReader(experiment_filename, bin_delimiter=args.bin_delimiter) as peptide_reader:
            psm_as_dict = []
            for psm in peptide_reader:
                psm_as_dict.append(psm)
        psm_data = pd.DataFrame(columns=data_store.DATAFRAME_COLUMNS)
        psm_data = psm_data.astype(data_store.DATAFRAME_DATATYPES)
        psm_data = pd.concat([psm_data, pd.DataFrame(psm_as_dict)], ignore_index=True)
        # psm_data.append(psm_as_dict, ignore_index=True)
        log(f'Parsed {len(psm_data.index)} PSMs from file.')
        spectrum_file_counter = 0

        ms_runs = psm_data["ms_run"].unique()
        for ms_run in ms_runs:
            try:
                ms_run_file = my_data_store.ms_run_address_book[ms_run]
            except KeyError:
                print(f'WARNING: No ms-run-file associated with ms run id "{ms_run}", skipping this file.')
                continue
            spectrum_file_counter += 1
            psm_data_of_current_ms_run = psm_data.loc[lambda df: df['ms_run'] == ms_run, :]

            # (2) sometimes, two psm's share the same ms2 spectrum. These are flagged here.
            flag_ambiguous_psms(psm_data_of_current_ms_run)

            # (3) read the ms_run file to grab currently unknown precursor ids for the known ms2 spectra as well as
            # collect all the ms1 spectra for subsampling, then: subsample.
            log(f'Now parsing ({spectrum_file_counter}/{len(ms_runs)}) spectrum files to find '
                f'{len(psm_data_of_current_ms_run.index)} MS1 precursors, "{os.path.basename(ms_run_file)}"...')
            ms1_spectra_hash = {}
            ms1_id_list = []
            ms1_spectra_list = []
            ms2_spectra_precursors = {}
            with pymzml.run.Reader(ms_run_file) as mzml_reader:
                for spectrum in mzml_reader:
                    if spectrum.ms_level == 1:
                        ms1_spectra_hash[spectrum.ID] = len(ms1_spectra_list)
                        ms1_id_list.append(spectrum.ID)
                        ms1_spectra_list.append(spectrum.peaks('centroided'))
                    elif spectrum.ms_level == 2:
                        ms2_spectra_precursors[spectrum.ID] = spectrum.selected_precursors
            psm_with_ms1_spectrum_count = 0
            for psm_index in psm_data_of_current_ms_run.index:
                try:
                    psm_with_ms1_spectrum_count += set_precursor(ms2_spectra_precursors[
                                                                 psm_data_of_current_ms_run.at[psm_index, 'psm_id']],
                                                             psm_data_of_current_ms_run, psm_index)
                except KeyError:
                    pass
            log(f'({psm_with_ms1_spectrum_count}/{len(psm_data_of_current_ms_run.index)}) PSMs have precursor MS1 spectrum')
            log(f'Now subsampling MS1 spectra for {len(psm_data_of_current_ms_run.index)} PSM patterns from '
                f'{len(ms1_spectra_list)} MS1 spectra, using {args.threads} threads for isotope analysis, '
                f'first prepping tasks for parallel execution...')

            pattern_data = pd.DataFrame(columns=data_store.DATAFRAME_COLUMNS)
            pattern_data = pattern_data.astype(data_store.DATAFRAME_DATATYPES)
            with ProcessPoolExecutor(max_workers=args.threads) as executor:
                peptides = psm_data_of_current_ms_run['peptide'].unique()
                tasks = [{'peptide': peptide,
                          'ms1_spectra_hash': ms1_spectra_hash,
                          'ms1_id_list': ms1_id_list,
                          'ms1_spectra_list': ms1_spectra_list,
                          'psms': psm_data_of_current_ms_run.loc[lambda df: df['peptide'] == peptide, :],
                          # 'psms': psm_data_of_current_ms_run[psm_data_of_current_ms_run['peptide'] == peptide],
                          } for peptide in peptides]
                patterns_reassigned_count = 0
                all_patterns = []
                for patterns in list(tqdm(executor.map(subsample_ms1_spectra, tasks, chunksize=25), total=len(tasks))):
                    for p in patterns:
                        if p['reassigned']:
                            patterns_reassigned_count += 1
                        del p['reassigned']
                    all_patterns.extend(patterns)
            pattern_data = pd.concat([pattern_data, pd.DataFrame(all_patterns)], ignore_index=True)
            # pattern_data = pattern_data.append(all_patterns, ignore_index=True)
            log(f'Subsampled and analyzed {len(pattern_data.index)} isotopic patterns, '
                f'reassigned {patterns_reassigned_count / len(pattern_data.index) * 100:.1f}% of patterns to '
                f'their nearest PSM. Now flagging overlap between patterns...')

            # (4) check for contamination (the same ms1 peaks shared by multiple peptide spectra
            with ProcessPoolExecutor(max_workers=args.threads) as executor:
                precursor_ids = pattern_data['pattern_precursor_id'].unique()
                tasks = [{'patterns from same spectrum': pattern_data[pattern_data['pattern_precursor_id'] == i]
                          } for i in precursor_ids]
                for indexes in list(tqdm(executor.map(
                        submit_for_overlap_detection, tasks, chunksize=25), total=len(tasks))):
                    for index in indexes:
                        pattern_data.at[index, 'flag_pattern_is_contaminated'] = True
            co = len(pattern_data.loc[lambda df: df["flag_pattern_is_contaminated"], :]) / len(pattern_data) * 100
            wo = len(pattern_data.loc[lambda df: df["flag_pattern_is_wobbly"], :]) / len(pattern_data) * 100
            log(f'Contamination found for {co:.1f}% of subsampled patterns. Wobbly patterns make up {wo:.1f}%.')
            ms_run_results_file = os.path.join(args.output_file,
                                               os.path.splitext(os.path.basename(ms_run_file))[0] + '.feather')
            log(f'Saving pandas dataframe with {len(pattern_data.index)} analyzed patterns to f{ms_run_results_file} '
                f'in feather format...')
            pattern_data.to_feather(ms_run_results_file)
    log(f'Done. Thanks for using calisp!')


if __name__ == "__main__":
    main()
