from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from pandas.core.interchange.dataframe_protocol import DataFrame

from calisp.log import log, VERSION

FILES = []
FFT_MAX_ERROR = 0.001
CRAP = False
FLAGS = ('flag_psm_has_low_confidence flag_psm_is_ambiguous flag_pattern_is_contaminated flag_pattern_is_wobbly '
         'flag_peptide_assigned_to_multiple_proteins flag_peptide_assigned_to_multiple_bins flag_peptide_mass_and_elements_undefined').split()

def parse_arguments():
    global FFT_MAX_ERROR, FLAGS, FILES, CRAP
    parser = ArgumentParser(description='calisp_filter_patterns.py. (C) Marc Strous, 2025')
    parser.add_argument('--result_file', required=True,  # type=argparse.FileType('r'),
                        help='[.feather] file or folder with [.feather] calisp output file(s).')
    parser.add_argument('--SIF', default=False, action='store_true', help='SIF (fingerprinting, natural abundances).')
    parser.add_argument('--SIP', default=False, action='store_true', help='SIP (probing, isotopes were added to the experiment).')
    parser.add_argument('--CRAP', default=False, action='store_true', help='Remove patterns of proteins of the CRAP database.')
    parser.add_argument('--flags', help='For expert use only: comma separated list of flags. Patterns that have one or more of '
                                                     'these flags will be filtered out.')
    parser.add_argument('--max_fft_error', help='For expert use only: Only keep patterns that fit the model better '
                                                             'than this value.')

    args = parser.parse_args()
    if args.SIF:
        FLAGS = 'flag_peptide_assigned_to_multiple_bins flag_peptide_mass_and_elements_undefined flag_pattern_is_contaminated'.split()
        log('Filters set for stable isotope fingerprinting:')
        log(f'  FFT maximum error = {FFT_MAX_ERROR}')
        log(f'  Flags: {", ".join(FLAGS)}')
    elif args.SIP:
        FFT_MAX_ERROR = 0.0
        log('Filters set for stable isotope probing:')
        log(f'  FFT maximum error = (not used)')
        log(f'  Flags: {", ".join(FLAGS)}')
    else:
        raise Exception('Need to specify --SIP or --SIF.')
    if args.CRAP:
        CRAP = True
    if args.max_fft_error:
        FFT_MAX_ERROR = float(args.max_fft_error)
        log(f'  (Expert override of default) FFT maximum error = {FFT_MAX_ERROR}')
    if args.flags:
        FLAGS = args.flags.split(',')
        log(f'  (Expert override of default) Flags: {", ".join(FLAGS)}')
    input_file = Path(args.result_file)
    if not input_file.exists():
        raise Exception(f'Input file {input_file} not found.')
    if input_file.is_file():
        output_file = Path(input_file.parent, input_file.stem + '.filtered.csv')
        FILES.append({'in': input_file, 'out': output_file})
    elif input_file.is_dir():
        output_dir = Path(input_file.parent, input_file.name + '.filtered')
        output_dir.mkdir(exist_ok=True)
        for f in input_file.glob('*.feather'):
            FILES.append({'in': f, 'out': Path(output_dir, f.stem + '.filtered.csv')})
    log('Files: ')
    for f in FILES:
        log(f'  {f["in"].parent.name}/{f["in"].name} => {f["out"].parent.name}/{f["out"].name}')


def filter_calisp_data(data:DataFrame, flags:list, fft_max_error:float, crap:bool):
    initial_pattern_count = len(data.index)
    if 'flag_peptide_mass_and_elements_undefined' in flags:
        data = data[data.flag_peptide_mass_and_elements_undefined != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_peptide_mass_and_elements_undefined"...')
    if 'flag_peptide_assigned_to_multiple_bins' in flags:
        data = data[data.flag_peptide_assigned_to_multiple_bins != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_peptide_assigned_to_multiple_bins"...')
    if 'flag_peptide_assigned_to_multiple_proteins' in flags:
        data = data[data.flag_peptide_assigned_to_multiple_proteins != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_peptide_assigned_to_multiple_proteins"...')
    if 'flag_psm_has_low_confidence' in flags:
        data = data[data.flag_psm_has_low_confidence != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_psm_has_low_confidence"...')
    if 'flag_psm_is_ambiguous' in flags:
        data = data[data.flag_psm_is_ambiguous != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_psm_is_ambiguous"...')
    if 'flag_pattern_is_contaminated' in flags:
        data = data[data.flag_pattern_is_contaminated != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_pattern_is_contaminated"...')
    if 'flag_pattern_is_wobbly' in flags:
        data = data[data.flag_pattern_is_wobbly != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_pattern_is_wobbly"...')
    if 'flag_peak_at_minus_one_pos' in flags:
        data = data[data.flag_peak_at_minus_one_pos != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_peak_at_minus_one_pos"...')
    if 'flag_peptide_contains_sulfur' in flags:
        data = data[data.flag_peptide_contains_sulfur != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_peptide_contains_sulfur"...')
    if 'flag_peptide_has_modifications' in flags:
        data = data[data.flag_peptide_has_modifications != True]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding "flag_peptide_has_modifications"...')
    if crap:
        data = data[~data.proteins.str.contains('^CRAP')]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding protein names containing "^CRAP"...')
    if fft_max_error > 0:
        data = data[data.error_fft < fft_max_error]
        log(f'{len(data.index)} ({len(data.index) / initial_pattern_count:.1%}) remaining after excluding patterns with a FFT fit worse than {fft_max_error}...')
    return data


def main():
    log(f'This is calisp_filter_patterns.py, version {VERSION}')
    args = parse_arguments()
    for f in FILES:
        data = pd.read_feather(f['in'])
        log(f'{f["in"].name}: Loaded {len(data.index)} patterns.')
        data = filter_calisp_data(data, flags=FLAGS, fft_max_error=FFT_MAX_ERROR, crap=CRAP)
        log(f'Writing filtered patterns to {f["out"]}...')
        log('=========')
        data.to_csv(f['out'])
        log(f'Thank you for using calisp_filter_patterns.py, version {VERSION}')

if __name__ == "__main__":
    main()
