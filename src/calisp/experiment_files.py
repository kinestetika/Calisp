import os
import re
import numpy as np
from calisp import element_count_and_mass_utils
from calisp import isotopic_pattern_utils
from calisp.log import log


UNBINNED = 'unbinned'
PATTERN_BOUNDARY_AA = re.compile('\[[A-Z]]')

def get_column_indexes(header_line, platform_target_columns):
    header_words = header_line.replace('"', '').strip().split('\t')
    return {c: header_words.index(c) for c in platform_target_columns}


def parse_psm_info_fragpipe(line, target_columns, unknown_modifications):
    return {
        'sequence': '',
        'above_threshold': False,
        'ambiguous': False,
        'modifications': [],
        'ms2_id': 0,
        'spectrum_file': '',
        'charge': 0,
        'm/z': 0.0,
        'proteins': {},
    }


def parse_psm_info_proteome_discoverer(line, target_columns, unknown_modifications):
    words = line[1:-1].split('"\t"')
    # strip mods from aminoacid sequence:
    peptide_aminoacid_sequence = words[target_columns['Annotated Sequence']].upper()
    peptide_aminoacid_sequence = re.sub(PATTERN_BOUNDARY_AA, '', peptide_aminoacid_sequence).strip('.')
    # strip extension from spectrum file
    spectrum_file_name = words[target_columns['Spectrum File']]
    (spectrum_file_base, extension) = os.path.splitext(spectrum_file_name)
    # parse modifications
    peptide_modifications = []
    peptide_modification_positions = []
    for modification_definition in words[target_columns['Modifications']].split(';'):
        modification_definition = modification_definition.strip()
        if not modification_definition:
            continue
        m = re.match('[A-Z]([0-9]+)[(](.+?)[)]', modification_definition)
        if not m:
            m = re.match('([CN])-Term(?:[(]Prot[)])*[(](.+?)[)]', modification_definition)
        if not m:
            peptide_modifications.append(-1)
            peptide_modification_positions.append(0)
            unknown_modifications.add(modification_definition)
            continue
        i = element_count_and_mass_utils.unimod_peptide_modifications_id(m.group(2))
        if i >= 0:
            pos = 0
            try:
                pos = int(m.group(1)) - 1
            except ValueError:
                if 'C' == m.group(1):
                    pos = len(peptide_aminoacid_sequence) - 1
            else:
                peptide_modifications.append(i)
                peptide_modification_positions.append(pos)
        else:
            peptide_modifications.append(-1)
            peptide_modification_positions.append(0)
            unknown_modifications.add(modification_definition)
    return {
        'sequence': peptide_aminoacid_sequence,
        'above_threshold': words[target_columns['Confidence']] == 'High',
        'ambiguous': 'unambiguous' != words[target_columns['PSM Ambiguity']].lower(),
        'modifications': peptide_modifications,
        'modification_positions': peptide_modification_positions,
        'ms2_id': int(words[target_columns['First Scan']]),
        'spectrum_file': spectrum_file_base,
        'charge': int(words[target_columns['Charge']]),
        'm/z': float(words[target_columns['m/z [Da]']]),
        'proteins': {p.strip() for p in words[target_columns['Protein Accessions']].replace("\t", " ").split(';')},
    }


PLATFORMS = {'proteome_discoverer': {'column_names': ['Annotated Sequence', 'Confidence', 'PSM Ambiguity', 'Modifications', 'First Scan',
                                                      'Spectrum File', 'Charge', 'm/z [Da]', 'Protein Accessions'],
                                     'parser': parse_psm_info_proteome_discoverer},
             'fragpipe':            {'column_names': ['Peptide', 'Expectation', 'Assigned Modifications', 'Spectrum', 'Charge',
                                                      'Calibrated Observed M/Z', 'Protein', 'Mapped Proteins'],
                                     'parser': parse_psm_info_fragpipe},
            }


class TabularExperimentFileReader:
    def __init__(self, experiment_file_name, bin_delimiter='_'):
        self.experiment_file_name = experiment_file_name
        self.bin_delimiter = bin_delimiter
        self.parse_bins = bin_delimiter != 'none'
        self.unknown_modifications = set()

    def __enter__(self):
        self.file = open(self.experiment_file_name)
        line = self.file.readline()
        for platform, platform_def in PLATFORMS.items():
            try:
                self.target_columns = get_column_indexes(line, platform_def['column_names'])
                self.parser = platform_def['parser']
                log(f'Experiment file "{self.experiment_file_name}" generated with {platform}')
                break
            except ValueError:
                pass
        if not self.target_columns:
            log(f'Experiment file "{self.experiment_file_name}" has unknown format. Aborting...')
            exit(1)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        return self

    def close(self):
        self.file.close()

    def __next__(self) -> {}:
        for line in self.file:
            psm_info = self.parser(line, self.target_columns, self.unknown_modifications)
            # parse bin names
            proteins = set()
            bins = set()
            if self.parse_bins:
                for protein in psm_info['proteins']:
                    bin = UNBINNED
                    if self.bin_delimiter in protein:
                        (bin, protein) = protein.split(self.bin_delimiter, 1)
                    proteins.add(protein)
                    bins.add(bin)
            # calculate peptide element composition, mass, and neutrons
            peptide_element_counts = np.zeros(5, dtype=np.int16)  # C, N, O, H, S
            try:
                peptide_mass = element_count_and_mass_utils.compute_peptide_mass(psm_info['sequence'],
                                                                                 psm_info['modifications'],
                                                                                 peptide_element_counts)
                psm_neutrons_in_ms2_spectrum = sum(
                    element_count_and_mass_utils.unimod_peptide_modifications_extra_neutrons(ptm)
                    for ptm in psm_info['modifications'])
                psm_neutrons_in_ms2_spectrum /= peptide_element_counts[isotopic_pattern_utils.ELEMENT_ROW_INDEX]
                peptide_mass_undefined = False
            except KeyError:
                peptide_element_counts = np.zeros(5, dtype=np.int16)
                peptide_mass_undefined = True
                psm_neutrons_in_ms2_spectrum = 0
                peptide_mass = 0

            return {'experiment': os.path.basename(self.experiment_file_name),
                    'ms_run': psm_info['spectrum_file'],
                    'bins': ' '.join(bins),
                    'proteins': ' '.join(proteins),
                    'peptide': f'{psm_info['sequence']} '
                               f'{",".join([element_count_and_mass_utils.unimod_list()[m] for m in psm_info["modifications"]])} '
                               f'{",".join(str(p) for p in psm_info["modification_positions"])}',
                    'psm_id': psm_info['ms2_id'],  # this refers to the id of the associated ms2 spectrum
                    'psm_mz': psm_info['m/z'],
                    'psm_charge': psm_info['charge'],
                    'psm_neutrons': psm_neutrons_in_ms2_spectrum,
                    'psm_rank': 0,
                    'psm_precursor_id': 0,
                    'psm_precursor_mz': 0,
                    'C': peptide_element_counts[0],
                    'N': peptide_element_counts[1],
                    'O': peptide_element_counts[2],
                    'H': peptide_element_counts[3],
                    'S': peptide_element_counts[4],
                    'peptide_mass': peptide_mass,
                    'flag_peptide_contains_sulfur': element_count_and_mass_utils.contains_sulfur(psm_info['sequence']),
                    'flag_peptide_has_modifications': len(psm_info['modifications']) > 0,
                    'flag_peptide_assigned_to_multiple_bins': len(bins) > 1,
                    'flag_peptide_assigned_to_multiple_proteins': len(proteins) > 1,
                    'flag_peptide_mass_and_elements_undefined': peptide_mass_undefined,
                    'flag_psm_has_low_confidence': not psm_info['above_threshold'],  # proteome discoverer only
                    'flag_psm_is_ambiguous': psm_info['ambiguous']}  # proteome discoverer only

        # end of file reached
        # log problems
        if len(self.unknown_modifications):
            print(f'Warning: Encountered {len(self.unknown_modifications)} unknown modifications: ' +
                  ' '.join(self.unknown_modifications))
        raise StopIteration
