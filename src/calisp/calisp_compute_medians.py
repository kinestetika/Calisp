from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from calisp.log import log, VERSION
from calisp import isotopic_pattern_utils

FILES = []
VOCABULARY = []
DELTA = False
RESULT_FILE = Path()
TARGET_ISOTOPE = 'C13'
BASE_ISOTOPE = 'C12'
FOR_EACH_PROTEIN = False

def parse_arguments():
    global FILES, VOCABULARY, DELTA, RESULT_FILE, TARGET_ISOTOPE, BASE_ISOTOPE, FOR_EACH_PROTEIN
    parser = ArgumentParser(description='calisp_compute_medians.py. (C) Marc Strous, 2025')
    parser.add_argument('--result_file', required=True,  # type=argparse.FileType('r'),
                        help='[.csv] file or folder with [.csv] calisp filtered output file(s).')
    parser.add_argument('--SIF', default=False, action='store_true', help='SIF (fingerprinting, natural abundances).')
    parser.add_argument('--SIP', default=False, action='store_true', help='SIP (probing, isotopes were added to the experiment).')
    parser.add_argument('--protein', default=False, action='store_true', help='Collect stats per protein instead of bin.')
    parser.add_argument('--vocabulary_file', default='', help='a tab delimited file defining any changes that need to be made to '
                                                                           '"proteins", "bins", "experiment" and "ms_run" (see README.md for examples)')
    parser.add_argument('--isotope', default='13C', choices=['13C', '14C', '15N', '17O', '18O', '2H', '3H', '33S',
                                                             '34S', '36S'],
                        help='The target isotope. Default: 13C. Only needed for computation of delta value (SIF)')
    parser.add_argument('--isotope_abundance_matrix', default=None,
                        help='To use a custom isotope abundance matrix. Only used for computation of delta value (SIF)')

    args = parser.parse_args()
    input_file = Path(args.result_file)
    if not input_file.exists():
        raise Exception(f'Input file {input_file} not found.')
    if input_file.is_file():
        FILES.append(input_file)
        RESULT_FILE = Path(input_file.parent, input_file.stem + '.stats.csv')
    elif input_file.is_dir():
        for f in input_file.glob('*.filtered.csv'):
            FILES.append(f)
        RESULT_FILE = Path(input_file, 'stats.csv')
    log('Files: ')
    for f in FILES:
        log(f'  {f}')
    if args.SIF:
        DELTA = True
    elif args.SIP:
        DELTA = False
    else:
        log('Need to specify --SIP or --SIF.')
        exit(1)
    TARGET_ISOTOPE = args.isotope
    BASE_ISOTOPE = {'13C': '12C', '14C': '12C', '15N': '14N', '17O': '16O', '18O': '16O',
                    '33S': '32S', '34S': '32S', '36S': '32S'}[args.isotope]
    (isotopic_pattern_utils.ELEMENT_ROW_INDEX, isotopic_pattern_utils.ISOTOPE_COLUMN_INDEX) = \
        {'13C': (0, 1), '14C': (0, 2),
         '15N': (1, 1),
         '17O': (2, 1), '18O': (2, 2),
         '2H': (3, 1), '3H': (3, 2),
         '33S': (4, 1), '34S': (4, 2), '36S': (4, 4)}[args.isotope]
    if args.isotope_abundance_matrix:
        isotopic_pattern_utils.DEFAULT_MATRIX_FILE = Path(args.isotope_abundance_matrix)
    isotopic_pattern_utils.ISOTOPE_MATRIX = isotopic_pattern_utils.load_isotope_matrix(None)
    isotopic_pattern_utils.NATURAL_ABUNDANCES = isotopic_pattern_utils.load_isotope_matrix(None)

    if args.vocabulary_file:
        vocabulary_file = Path(args.vocabulary_file)
        if not input_file.exists():
            raise Exception(f'Input file {input_file} not found.')
        with open(vocabulary_file) as reader:
            for line in reader:
                try:
                    words = line.split('\t')
                    words = {'src_key': words[0],
                             'src_value': words[1],
                             'target_key': words[2],
                             'target_value': words[3]}
                    VOCABULARY.append(words)
                except:
                    log('The vocabulary file should contain four names per row, separated by a tab.')
                    log(f'Offending line: "{line}"')
                    log('Please fix the vocabulary and run again.')
                    exit(1)
    FOR_EACH_PROTEIN = args.protein


def parse_data(filenames):
    data = []
    for f in filenames:
        new_data = pd.read_csv(f)
        log(f'{f.name}: Parsed {len(new_data.index)} patterns.')
        data.append(new_data)
    data = pd.concat(data)
    log(f'Parsed {len(data.index)} patterns in total.')
    return data


def apply_vocabulary(data, vocabulary):
    for entry in vocabulary:
        affected_data = data[entry["src_key"]] == entry["src_value"]
        data.loc[affected_data, entry["target_key"]] = entry["target_value"]
    return data


def count_unique_proteins(data):
    proteins = list(data['proteins'].unique())
    proteins = set(' '.join(proteins).split())
    return len(proteins)


def get_topic_stats(for_each_protein, topic, experiment, ms_run, target_column, target_column_name, data):
    if for_each_protein:
        return {'protein(s)': topic,
                'bins(s)': ' '.join(set(' '.join(data.bins.unique()).split())),
                'experiment': experiment,
                'ms_run': ms_run,
                'summed intensity': data['pattern_total_intensity'].sum(),
                '# unique peptides': len(data['peptide'].unique()),
                '# PSMs': len(data['psm_id'].unique()),
                '# patterns': len(data.index),
                '#spectra': len(data.index),
                target_column_name: data[target_column].median(),
                'lower_quantile': data[target_column].quantile(0.25),
                'upper_quantile': data[target_column].quantile(0.75),
                'mean_psm_neutron_count': float(data['psm_neutrons'].sum()) / len(data.index)}
    else:
        return {'bin': topic,
                'experiment': experiment,
                'ms_run': ms_run,
                'summed intensity': data['pattern_total_intensity'].sum(),
                '# unique proteins': count_unique_proteins(data),
                '# unique peptides': len(data['peptide'].unique()),
                '# PSMs': len(data['psm_id'].unique()),
                '# patterns': len(data.index),
                '#spectra': len(data.index),
                target_column_name: data[target_column].median(),
                'lower_quantile': data[target_column].quantile(0.25),
                'upper_quantile': data[target_column].quantile(0.75),
                'mean_psm_neutron_count': float(data['psm_neutrons'].sum()) / len(data.index)}


def aggregate_stats_for_topic(for_each_protein, topic, target_column, target_column_name, all_stats, new_data):
    for experiment in new_data.experiment.unique():
        experiment_data = new_data[new_data.experiment == experiment]
        for ms_run in experiment_data.ms_run.unique():
            ms_run_data = experiment_data[experiment_data.ms_run == ms_run]
            ms_run_stats = get_topic_stats(for_each_protein, topic, experiment, ms_run, target_column,
                                           target_column_name, ms_run_data)
            if all_stats is not None:
                all_stats = pd.concat([all_stats, pd.DataFrame([ms_run_stats])], ignore_index=True)
            else:
                all_stats = pd.DataFrame([ms_run_stats])
    return all_stats


def compute_medians(data, delta, for_each_protein):
    if delta:
        target_column = "ratio_fft"
        target_column_name = f"median delta{TARGET_ISOTOPE} (per mille)"
        #data["ratio_fft"] = (data["ratio_fft"] / (0.011056585166521 / 0.988943414833479) - 1) * 1000
        data["ratio_fft"] = (data["ratio_fft"] / (isotopic_pattern_utils.ISOTOPE_MATRIX[isotopic_pattern_utils.ELEMENT_ROW_INDEX]
                                                                                       [isotopic_pattern_utils.ISOTOPE_COLUMN_INDEX] /
                                                  isotopic_pattern_utils.ISOTOPE_MATRIX[isotopic_pattern_utils.ELEMENT_ROW_INDEX][0]) - 1) * 1000
    else:
        target_column = "ratio_na"
        target_column_name = f"median ratio {TARGET_ISOTOPE}/{BASE_ISOTOPE} (-)"

    all_stats = None
    if for_each_protein:
        for protein in data.proteins.unique():
            protein_data = data[data.proteins == protein]
            all_stats = aggregate_stats_for_topic(for_each_protein, protein, target_column, target_column_name, all_stats, protein_data)
    else:
        for bin in data.bins.unique():
            bin_data = data[data.bins == bin]
            all_stats = aggregate_stats_for_topic(for_each_protein, bin, target_column, target_column_name, all_stats, bin_data)
    return all_stats


def main():
    log(f'This is calisp_compute_medians.py, version {VERSION}')
    parse_arguments()
    data = parse_data(FILES)
    data = apply_vocabulary(data, VOCABULARY)
    stats = compute_medians(data, DELTA, FOR_EACH_PROTEIN)
    stats.to_csv(RESULT_FILE, index = False)


if __name__ == "__main__":
    main()
