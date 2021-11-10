import os
import re
import numpy as np
from src import element_count_and_mass_utils as util, spectrum_analysis_utils

UNBINNED = 'unbinned'


class PeptideSpectrumMatchFileReader:
    def __init__(self, experiment_file_name, bin_delimiter ='_'):
        self.experiment_file_name = experiment_file_name
        self.bin_delimiter = bin_delimiter
        self.unknown_modifications = set()

    def __enter__(self):
        try:
            self.file = open(self.experiment_file_name)
            line = self.file.readline()
            line = line[1:-1]
            words = line.split("\"\t\"")
            self.columns = {'Annotated Sequence': -1,
                            'Confidence': -1,
                            'PSM Ambiguity': -1,
                            'Modifications': -1,
                            'First Scan': -1,
                            'Spectrum File': -1,
                            'Charge': -1,
                            'm/z [Da]': -1,
                            'Protein Accessions': -1
                           }
            for k in self.columns.keys():
                self.columns[k] = words.index(k)
            return self
        except ValueError:
            raise Exception(f"Missing column in spectrum match file {self.experiment_file_name}")

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        return self

    def close(self):
        self.file.close()

    def __next__(self) -> {}:
        for line in self.file:
            # line = self.file.readline()
            line = line[1:-1]
            words = line.split("\"\t\"")
            psm_spectrum_id = int(words[self.columns['First Scan']])
            psm_is_above_threshold = words[self.columns['Confidence']] == 'High'
            psm_charge = int(words[self.columns['Charge']])
            psm_mz = float(words[self.columns['m/z [Da]']])
            psm_neutrons_in_ms2_spectrum = 0
            psm_is_ambiguous = 'unambiguous' != words[self.columns['PSM Ambiguity']].lower()

            peptide_aminoacid_sequence = words[self.columns['Annotated Sequence']].upper()
            peptide_modifications = []
            peptide_modification_positions = []
            ms_run_name = words[self.columns['Spectrum File']]
            (ms_run_name, extension) = os.path.splitext(ms_run_name)
            peptide_protein_ids = set()
            peptide_bin_ids = set()
            peptide_protein_names = set()
            peptide_bin_names=set()
            for protein_acc in words[self.columns['Protein Accessions']].replace("\t", " ").split(';'):
                protein_acc = protein_acc.strip()
                bin_name = ''
                protein_name = ''
                if self.bin_delimiter in protein_acc:
                    (bin_name, protein_name) = protein_acc.split(self.bin_delimiter, 1)
                else:
                    protein_name = protein_acc
                    bin_name = UNBINNED
                #peptide_bin_ids.add(self.data_store.add_to_store(bin_name, self.data_store.bin_hash,
                #                                                 self.data_store.bins))
                #peptide_protein_ids.add(self.data_store.add_to_store(protein_name, self.data_store.protein_hash,
                #                                                 self.data_store.proteins))
                peptide_protein_names.add(protein_name)
                peptide_bin_names.add(bin_name)

            peptide_mass_undefined = False
            for modification_definition in words[self.columns['Modifications']].split(';'):
                modification_definition = modification_definition.strip()
                if not modification_definition:
                    continue
                m = re.match('[A-Z]([0-9]+)[(](.+?)[)]', modification_definition)
                if not m:
                    m = re.match('([CN])-Term(?:[(]Prot[)])*[(](.+?)[)]', modification_definition)
                if not m:
                    peptide_mass_undefined = True
                    self.unknown_modifications.add(modification_definition)
                    continue
                i = util.unimod_peptide_modifications_id(m.group(2))
                if i >= 0:
                    pos = 0
                    try:
                        pos = int(m.group(1))-1
                    except ValueError:
                        if 'C' == m.group(1):
                            pos = len(peptide_aminoacid_sequence) - 1
                    neutrons = util.unimod_peptide_modifications_extra_neutrons(i)
                    if neutrons > 0:
                        psm_neutrons_in_ms2_spectrum += util.unimod_peptide_modifications_extra_neutrons(i)
                    else:
                        peptide_modifications.append(i)
                        peptide_modification_positions.append(pos)
                else:
                    peptide_mass_undefined = True
                    self.unknown_modifications.add(modification_definition)

            peptide_element_counts = np.zeros(5, dtype=np.int16)  # C, N, O, H, S
            try:
                peptide_mass = util.compute_peptide_mass(peptide_aminoacid_sequence,
                                                       peptide_modifications,
                                                       peptide_element_counts)
                psm_neutrons_in_ms2_spectrum /= peptide_element_counts[spectrum_analysis_utils.ELEMENT_ROW_INDEX]
            except KeyError:
                peptide_element_counts = np.zeros(5, dtype=np.int16)
                peptide_mass = 0
                peptide_mass_undefined = True
                psm_neutrons_in_ms2_spectrum = 0
            except ValueError:
                print(peptide_aminoacid_sequence)
                exit(1)

            return {'experiment': os.path.basename(self.experiment_file_name),
                    'ms_run': ms_run_name,
                    'bins': ' '.join(peptide_bin_names),
                    'proteins': ' '.join(peptide_protein_names),
                    'peptide': f'{peptide_aminoacid_sequence} '
                               f'{[util.unimod_list()[m] for m in peptide_modifications]} '
                               f'{peptide_modification_positions}',
                    'psm_id': psm_spectrum_id,
                    'psm_mz': psm_mz,
                    'psm_charge': psm_charge,
                    'psm_neutrons': psm_neutrons_in_ms2_spectrum,
                    'psm_rank': 0,
                    'psm_precursor_id': 0, # precursor_ms1_spectrum_id
                    'psm_precursor_mz': 0., # precursor_ms1_spectrum_mz
                    'C': peptide_element_counts[0], #C'
                    'N': peptide_element_counts[1], #N
                    'O': peptide_element_counts[2], #O
                    'H': peptide_element_counts[3], #H
                    'S': peptide_element_counts[4], #S
                    'peptide_mass': peptide_mass,
                    'flag_peptide_contains_sulfur': util.contains_sulfur(peptide_aminoacid_sequence),
                    'flag_peptide_has_modifications': len(peptide_modifications) > 0,
                    'flag_peptide_assigned_to_multiple_bins': len(peptide_bin_names) > 1,
                    'flag_peptide_assigned_to_multiple_proteins': len(peptide_protein_names) > 1,
                    'flag_peptide_mass_and_elements_undefined': peptide_mass_undefined,
                    'flag_psm_has_low_confidence': not psm_is_above_threshold,
                    'flag_psm_is_ambiguous': psm_is_ambiguous}

        # end of file reached
        # log problems
        if len(self.unknown_modifications):
            print(f'Warning: Encountered {len(self.unknown_modifications)} unknown modifications, '
                    f'such as {self.unknown_modifications}')
        raise StopIteration
