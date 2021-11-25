import os.path
import pymzml
from pymzml.run import Reader

DATAFRAME_DATATYPES = {'experiment': str,
                       'ms_run': str,
                       'bins': str,
                       'proteins': str,
                       'peptide': str,
                       'peptide_mass': float,
                       'C': int,
                       'N': int,
                       'O': int,
                       'H': int,
                       'S': int,
                       'psm_id': int,
                       'psm_mz': float,
                       'psm_charge': int,
                       'psm_neutrons': float,
                       'psm_rank': int,
                       'psm_precursor_id': int,
                       'psm_precursor_mz': float,
                       'spectrum_charge': int,
                       'spectrum_precursor_id': int,
                       'spectrum_total_intensity': float,
                       'spectrum_peak_count': int,
                       'spectrum_median_peak_spacing': float,
                       'spectrum_mass_irregularity': float,
                       'ratio_na': float,
                       'ratio_fft': float,
                       'error_fft': float,
                       'error_clumpy': float,
                       'flag_peptide_contains_sulfur': bool,
                       'flag_peptide_has_modifications': bool,
                       'flag_peptide_assigned_to_multiple_bins': bool,
                       'flag_peptide_assigned_to_multiple_proteins': bool,
                       'flag_peptide_mass_and_elements_undefined': bool,
                       'flag_psm_has_low_confidence': bool,
                       'flag_psm_is_ambiguous': bool,
                       'flag_spectrum_is_contaminated': bool,
                       'flag_spectrum_is_wobbly': bool,
                       'flag_peak_at_minus_one_pos': bool}

DATAFRAME_COLUMNS = ['experiment', 'ms_run', 'bins', 'proteins', 'peptide', 'peptide_mass', 'C', 'N', 'O', 'H', 'S',
                     'psm_id', 'psm_mz', 'psm_charge', 'psm_neutrons', 'psm_rank', 'psm_precursor_id',
                     'psm_precursor_mz', 'spectrum_charge', 'spectrum_precursor_id', 'spectrum_total_intensity',
                     'spectrum_peak_count', 'spectrum_median_peak_spacing', 'spectrum_mass_irregularity',
                     'ratio_na', 'ratio_fft', 'error_fft', 'error_clumpy',
                     'flag_peptide_contains_sulfur', 'flag_peptide_has_modifications',
                     'flag_peptide_assigned_to_multiple_bins', 'flag_peptide_assigned_to_multiple_proteins',
                     'flag_peptide_mass_and_elements_undefined', 'flag_psm_has_low_confidence', 'flag_psm_is_ambiguous',
                     'flag_spectrum_is_contaminated', 'flag_spectrum_is_wobbly', 'flag_peak_at_minus_one_pos']

for i in range(20):  # intensities of the spectrum's peaks
    column_name = f'i{i}'
    DATAFRAME_COLUMNS.append(column_name)
    DATAFRAME_DATATYPES[column_name] = float

for m in range(20):  # masses of the spectrum's peaks
    column_name = f'm{m}'
    DATAFRAME_COLUMNS.append(column_name)
    DATAFRAME_DATATYPES[column_name] = float

for c in range(1, 7):  # clumpiness predictions
    column_name = f'c{c}'
    DATAFRAME_COLUMNS.append(column_name)
    DATAFRAME_DATATYPES[column_name] = float


def _is_mzid_file(file):
    return False


def _is_mzml_file(file):
    (file_name, extension) = os.path.splitext(file)
    return extension.lower() == '.mzml'


def _get_mzml_file_runid(file):
    parser: Reader = pymzml.run.Reader(file)
    run_id = parser.info.get('run_id')
    parser.close()
    return run_id


def _is_target_spectrum_match_file(f):
    try:
        file = open(f)
        line = file.readline()
        outcome = 'Annotated Sequence' in line and 'PSM Ambiguity' in line and 'm/z [Da]' in line
        file.close()
        return outcome
    except FileNotFoundError:
        return False
    except UnicodeDecodeError:
        return False
    except IsADirectoryError:
        return False


class DataStore:
    def __init__(self):
        self.experiments = []
        self.ms_runs = []
        self.ms_run_address_book = {}
        self.bins = ['unbinned']
        self.proteins = []
        self.peptides = []
        self.peptide_spectrum_matches = []
        self.spectra = []
        # the following are only used to speed up store construction
        self.bin_hash = {}
        self.protein_hash = {}
        self.peptide_hash = {}

    def set_scope(self, peptide_file_stub, spectrum_file_stub):
        if os.path.isfile(peptide_file_stub):
            self.experiments.append(peptide_file_stub)
        elif os.path.isdir(peptide_file_stub):
            for f in os.listdir(peptide_file_stub):
                f = os.path.join(peptide_file_stub, f)
                if _is_target_spectrum_match_file(f) or _is_mzid_file(f):
                    self.experiments.append(f)
        else:
            raise FileNotFoundError(f'Could not read or find provided peptide file: {peptide_file_stub!r}')
        if os.path.isfile(spectrum_file_stub):
            self.ms_runs.append(spectrum_file_stub)
        elif os.path.isdir(spectrum_file_stub):
            for f in os.listdir(spectrum_file_stub):
                f = os.path.join(spectrum_file_stub, f)
                if _is_mzml_file(f):
                    self.ms_runs.append(f)
        for f in self.ms_runs:
            (file_name, extension) = os.path.splitext(f)
            self.ms_run_address_book[file_name] = f
            self.ms_run_address_book[_get_mzml_file_runid(f)] = f
            self.ms_run_address_book[f] = f

    def get_ms_run_file_index(self, ms_run_file):
        return self.ms_runs.index(self.ms_run_address_book[ms_run_file])

    def get_experiment_index(self, experiment_file_name):
        return self.experiments.index(experiment_file_name)
