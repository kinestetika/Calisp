import numpy as np
import re
from calisp import unimod

AMINOACID_COUNT: int = 20
AMINOACID_ELEMENT_COUNTS = np.array([  # C, N, O, H, S
                                    [2, 1, 2, 5, 0],  # G
                                    [5, 1, 2, 9, 0],  # P
                                    [3, 1, 2, 7, 0],  # A
                                    [5, 1, 2, 11, 0],  # V
                                    [6, 1, 2, 13, 0],  # L
                                    [6, 1, 2, 13, 0],  # I
                                    [5, 1, 2, 11, 1],  # M
                                    [3, 1, 2, 7, 1],  # C
                                    [9, 1, 2, 11, 0],  # F
                                    [9, 1, 3, 11, 0],  # Y
                                    [11, 2, 2, 12, 0],  # W
                                    [6, 3, 2, 9, 0],  # H
                                    [6, 2, 2, 14, 0],  # K
                                    [6, 4, 2, 14, 0],  # R
                                    [5, 2, 3, 10, 0],  # Q
                                    [4, 2, 3, 8, 0],  # N
                                    [5, 1, 4, 9, 0],  # E
                                    [4, 1, 4, 7, 0],  # D
                                    [3, 1, 3, 7, 0],  # S
                                    [4, 1, 3, 9, 0]], dtype=np.int64)   # T
WATER_ELEMENT_COUNTS = np.array([0, 0, 1, 2, 0], dtype=np.int64)
AMINOACID_IDS = 'G P A V L I M C F Y W H K R Q N E D S T'.split()
AMINOACID_INDEXES = {x: AMINOACID_IDS.index(x) for x in AMINOACID_IDS}
ELEMENT_MONOISOTOPIC_MASSES = np.array([12.0, 14.0030740048, 15.99491461956, 1.0078250322, 31.972070999],
                                       dtype=np.float64)
AMINOACID_NAMES = 'Glycine(G) Proline(P) Alanine(A) Valine(V) Leucine(L) Isoleucine(I) Methionine(M) Cysteine(C) ' \
                  'Phenylalanine(F) Tyrosine(Y) Tryptophan(W) Histidine(H) Lysine(K) Arginine(R) Glutamine(Q) ' \
                  'Asparagine(N) Glutamate(E) Aspartate(D) Serine(S) Threonine(T)'.split()
AMINOACID_INDEXED_NAMES = {x: AMINOACID_NAMES[AMINOACID_IDS.index(x)] for x in AMINOACID_IDS}
PROTON_MASS = 1.0072765
NEUTRON_MASS_SHIFT = 1.002

__ALL_UNIMOD_PEPTIDE_MODIFICATIONS_ELEMENT_COUNTS = []
__ALL_UNIMOD_PEPTIDE_MODIFICATIONS_EXTRA_NEUTRONS = []
__ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS = []


def load_unicode():
    for line in unimod.UNIMOD.splitlines():
        words = line.split("\t")
        if len(words) < 2:
            continue
        m = re.search('C([0-9+-]+)H([0-9+-]+)N([0-9+-]+)O([0-9+-]+)S([0-9+-]+)', words[0])
        if m:
            unimod_id = words[1].strip()
            element_counts = np.empty(5, dtype=np.int16)
            element_counts[0] = int(m.group(1))
            element_counts[1] = int(m.group(3))
            element_counts[2] = int(m.group(4))
            element_counts[3] = int(m.group(2))
            element_counts[4] = int(m.group(5))
            extra_neutrons = np.int16(0)
            __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_ELEMENT_COUNTS.append(element_counts)
            __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS.append(unimod_id)
            if len(words) > 2:
                extra_neutrons = np.int16(int(words[2]))
            __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_EXTRA_NEUTRONS.append(extra_neutrons)


def unimod_list():
    if not __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS:
        load_unicode()
    return __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS


def unimod_peptide_modifications_element_counts(i):
    if not __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS:
        load_unicode()
    return __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_ELEMENT_COUNTS[i]


def unimod_peptide_modifications_extra_neutrons(i):
    if not __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS:
        load_unicode()
    return __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_EXTRA_NEUTRONS[i]


def unimod_peptide_modifications_id(unimod_id):
    if not __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS:
        load_unicode()
    for i in range(len(__ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS)):
        if __ALL_UNIMOD_PEPTIDE_MODIFICATIONS_IDS[i] == unimod_id:
            return i
    return -1


def compute_peptide_mass(aminoacid_sequence: str, peptide_modifications, element_counts=np.zeros(5, dtype=np.int16)):
    for i in range(5):
        element_counts[i] += WATER_ELEMENT_COUNTS[i]
    for aa in aminoacid_sequence.upper():
        aa_element_counts = AMINOACID_ELEMENT_COUNTS[AMINOACID_INDEXES[aa]]
        for i in range(5):
            element_counts[i] += aa_element_counts[i]
            element_counts[i] -= WATER_ELEMENT_COUNTS[i]
    for modification in peptide_modifications:
        modification_element_counts = unimod_peptide_modifications_element_counts(modification)
        if modification_element_counts is None:
            raise Exception(f"unknown modification {modification} in {aminoacid_sequence} "
                            f"- molecular mass computation failed")
        for i in range(5):
            element_counts[i] += modification_element_counts[i]
    mass = np.float64(0)
    for i in range(5):
        mass += element_counts[i] * ELEMENT_MONOISOTOPIC_MASSES[i]
    return mass


def get_element_counts_at_pos(aminoacid_sequence: str, peptide_modifications, peptide_modification_positions, pos):
    aa = aminoacid_sequence[pos]
    if aa not in AMINOACID_INDEXES:
        raise Exception(f"unknown aminoacid {aa} in {aminoacid_sequence} - elements at position computation failed")
    aa_element_counts = AMINOACID_ELEMENT_COUNTS[AMINOACID_INDEXES[aa]]
    element_counts = np.zeros(5, dtype=np.int64)
    for i in range(5):
        element_counts[i] = aa_element_counts[i]
    # eliminate water from peptide bond, except at C-terminus
    if pos < len(aminoacid_sequence)-1:
        for i in range(5):
            element_counts[i] -= WATER_ELEMENT_COUNTS[i]
    for j in range(len(peptide_modifications)):
        if peptide_modification_positions[j] == pos:
            modification_element_counts = unimod_peptide_modifications_element_counts(peptide_modifications[j])
            if not modification_element_counts:
                raise Exception(f"unknown modification {peptide_modifications[j]} in {aminoacid_sequence} "
                                f"- molecular mass computation at pos {pos} failed")
            for i in range(5):
                element_counts[i] += modification_element_counts[i]
    return element_counts


def contains_sulfur(aminoacid_sequence):
    return 'C' in aminoacid_sequence or 'M' in aminoacid_sequence
