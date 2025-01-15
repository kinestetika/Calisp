## Calisp.py, version 3.1

Calisp.py (Calgary approach to isotopes in proteomics) is a program that estimates isotopic composition (e.g. 13C/12C,
delta13C, 15N/14N etc) of peptides from proteomics mass spectrometry data. Input data consist of mzML files and 
files with peptide spectrum matches.

Calisp was originally developed in Java. This newer, python version is much more concise, it consists of only six .py 
files, ~1,000 lines of code, compared to about fifty files and ~10,000 lines of code for the Java program. In addition,
the Java program depends on another program, mcl, whereas calisp.py is purely python, which is easy to install, see below.
The conciseness of the code makes the python version more transparent, easier to maintain and easier to further develop. 
I consider calisp.py is the successor of the Java version. One of the reasons to shift to python is the possibility 
to more effectively develop machine learning approaches to filter out noisy isotopic patterns. Future work will explore
that possibility.

Benchmarking of Calisp.py has been completed. It works well, benchmarking procedures and outcomes are shared in the 
"benchmarking" folder. Parsing of .mzid files still needs to be implemented.

Calisp.py depends on numpy, scipy, pandas, tqdm, [pymzml](https://pymzml.readthedocs.io/en/latest/intro.html), pyarrow. 
These will be installed automatically by the pip command below. 
Calisp outputs the data as a Pandas DataFrame saved in (binary) [feather](https://arrow.apache.org/docs/python/feather.html) format.
Each row contains a single isotopic pattern, with column definitions listed below.
From there, the user can for example explore and visualize the results in a [Jupyter notebook](https://jupyter.org/). For that, a
tutorial wil be provided. In the meantime, you can check out the two notebooks I created for benchmarking calisp.py,
which are in the benchmarking folder.

Compared to previous versions of calisp, the workflow has been simplified. Calisp.py does not filter out any isotopic 
patterns, or adds up isotopic patterns to reduce noise - like the Java version does. It simply estimates the ratio for 
the target isotopes (e.g. 13C/12C) for every isotopic pattern it can subsample. It estimates this ratio based on neutron 
abundance and using fast fourier transforms. The former applies to stable isotope probing experiments. The latter applies 
to natural abundances, or to isotope probing experiments with very little added label (e.g. using substrates with <1% 
additional 13C). The motivation for omitting filtering is that keeping all subsampled isotopic patterns, including bad 
ones, will enable training of machine learning classifiers. Also, because it was shown that the median provides better 
estimates for species in microbial communities than the mean, adding up isotopic patterns to improve precision has lost 
its purpose. There is more power (and sensitivity) in numbers.

Because no data are filtered out and no isotopic patterns get added up, calisp.py, analyzes at least ten times as many 
isotopic patterns compared to the Java version. That means calisp.py is about ten times slower, it takes about 5-10 min 
per .mzML file on a Desktop computer. The user can perform filtering of the isotopic patterns using the Jupyter notebook 
as desired. For natural abundance data, it works well to only use those spectra that have a FFT fitting error 
("error_fft") of less than 0.001. Note that this threshold is less stringent th8en thew one used by the java program.

**Installation:**

>python -m virtualenv calisp

>source calisp/bin/activate

>pip install --upgrade calisp

If you would like to explore calisp results in Jupyter notebooks, run the following command instead:

>pip install --upgrade calisp jupyter matplotlib jinja2

**Usage:**

>source calisp/bin/activate

>calisp.py --spectrum_file [path to .mzML file or folder with .mzML files] --peptide_file [path to .peptideSpectrumMatch file 
 or folder with .PeptideSpectrumMatch files] --output_file [folder where calisp.py will save results files] --threads [# of 
 threads used, default 4] --isotope [15N, 3H etc, default 13C], --bin_delimiter [character that separates the bin ID from the
 remainder of protein IDs, default '_'], --mass_accuracy [accuracy of peak m/z identifications, default 10 ppm]
 --compute_clumps [use only if you want to compute clumpiness] --isotope_abundance_matrix [path to file with isotope matrix,
 a default file is included with calisp]

**Column names of the Pandas DataFrame created by calisp.py:**

In the saved dataframe, each row contains one isotopic pattern, defined by the following columns:
```
experiment             filename of the peptide spectrum match (psm) file
ms_run                 filename of the .mzml file
bins                   bin/mag ids, separated by commas. Calisp expects the protein ids in the psm 
                       file to consist of two parts, separated by a delimiter (_ by default). The 
                       first part is the bin/mag id, the second part the protein id
proteins               the ids of the proteins associated with the pattern (without the bin id)
peptide                the aminoacid sequence of the peptide
peptide_mass           the mass of the peptide
C                       # of carbon atoms in the peptide
N                       # of nitrogen atoms in the peptide
O                       # of oxygen atoms in the peptide
H                       # of hydrogen atoms in the peptide
S                       # of sulfur atoms in the peptide
psm_id                  psm id
psm_mz                  psm m over z
psm_charge              psm charge
psm_neutrons            number of neutrons inferred from custom 'neutron' modifications 
psm_rank                rank of the psm
psm_precursor_id        id of the ms1 spectrum that was the source of the psm 
psm_precursor_mz        mass over charge of the precursor of the psm
pattern_charge          charge of the pattern
pattern_precursor_id    id of the ms1 spectrum that was the source of the pattern
pattern_total_intensity total intensity of the pattern
pattern_peak_count      # of peaks in the pattern
pattern_median_peak_spacing medium mass difference between a pattern's peaks
spectrum_mass_irregularity  a measure for the standard deviation in the mass difference between a
                            pattern's peaks
ratio_na                the estimated isotope ratio inferred from neutron abundance (sip 
                        experiments) 
ratio_fft               the estimated isotope ratio inferred by the fft method (natural 
                        isotope abundances)
error_fft               the remaining error after fitting the pattern with fft
error_clumpy            the remaining error after fitting the pattern with the clumpy carbon 
                        method
flag_peptide_contains_sulfur true if peptide contains sulfur
flag_peptide_has_modifications true if peptide has no modifications
flag_peptide_assigned_to_multiple_bins true if peptide is associated with multiple proteins from 
                                       different bins/mags
flag_peptide_assigned_to_multiple_proteins true if peptide is associated with multiple proteins
flag_peptide_mass_and_elements_undefined true if peptide has unknown mass and elemental 
                                         composition
flag_psm_has_low_confidence true if psm was flagged as having low confidence (peptide identity 
                            uncertain)
flag_psm_is_ambiguous   true if psm could not be assigned with certainty
flag_pattern_is_contaminated true if multiple patterns have one or more shared peaks
flag_pattern_is_wobbly true if pattern_median_peak_spacing exceeds a treshold
flag_peak_at_minus_one_pos  true if a peak was detected immediately before the monoisotopic peak,
                            could indicate overlap with another pattern
i0 - i19                the intensities of the first 20 peaks of the pattern  
m0 - m19                the masses of the first 20 peaks of the pattern
c1 - c6                 contributions of clumps of 1-6 carbon to ratio_na. These are the 
                        outcomes of the clumpy carbon
                        model. These results are only meaningful if the biomass was labeled to 
                        saturation.
```
**Using a custom isotope abundance matrix:**
When estimating isotopic content for nitrogen, oxygen, hydrogen and sulfur, the estimates will be strongly
affected by isotope abundances of carbon. For example, often biomass is slightly depleted in 13C compared to
the inorganic reference (Vienna Pee Dee Belemnite). However, Calisp will assume the assumed 13C content to be
correct and will compensate the difference by reducing the content of the target isotope, which may result in
estimating a negative content for the target isotope (which is physically impossible). To overcome this issue,
you can provide a custom isotope matrix with the actual or estimated 13C content of your samples (or other changes
you may wish to make). The isotope matrix file should be formatted as follows, with elements on rows and isotopes
(+0, +1, +2, ... extra neutrons) on columns:
```
0.988943414833479 0.011056585166521 0.0         0.0 0.0    0.0 0.0 # C
0.996323567       0.003676433       0.0         0.0 0.0    0.0 0.0 # N
0.997574195       0.00038           0.002045805 0.0 0.0    0.0 0.0 # O
0.99988           0.00012           0.0         0.0 0.0    0.0 0.0 # H
0.9493            0.0076            0.0429      0.0 0.0002 0.0 0.0 # S
# VPDB standard 13C/12C = 0.0111802 in Isodat software
# see also https://www.webelements.com/sulfur/isotopes.html
# see also http://iupac.org/publications/pac/pdf/2003/pdf/7506x0683.pdf
```
**calisp_filter_patterns.py**
Use this script after running calisp.py to filter the detected patterns and to convert your data from binary .feather
format to .csv.

Example usage: calisp_filter_patterns.py --result_file [calisp result file or dir] --SIF

You need to specify --SIP or --SIF.

Arguments:
    --result_file       The .feather file generated by calisp, or a dir containing multiple .feather files.
    --SIF               Apply benchmarked filters for stable isotope fingerprinting data.
    --SIP               Apply benchmarked filters for stable isotope probing data.
    --CRAP              Remove proteins of the CRAP database.
    --flags             (Expert use:) Specify a custom list of flags for filtering patterns. Options are:
                            flag_psm_has_low_confidence*
                            flag_psm_is_ambiguous*
                            flag_pattern_is_contaminated*$
                            flag_pattern_is_wobbly*
                            flag_peptide_assigned_to_multiple_proteins*
                            flag_peptide_assigned_to_multiple_bins*$
                            flag_peptide_mass_and_elements_undefined*$
                            flag_peptide_has_modifications
                            flag_peptide_contains_sulfur
                        Flags with * are default when using --SIP. Flags with $ are default when using --SIF
    --max_fft_error     (Expert use:) Specify a custom maximum value for the allowable fft error. This only makes sense
                        for SIF data. The benchmarked default when using --SIF is 0.001.


Calisp was developed using [PyCharm community edition](https://www.jetbrains.com/pycharm/).

**Please cite:**

Kleiner M, Dong X, Hinzke T, Wippler J, Thorson E, Mayer B, Strous M (2018) A metaproteomics method to determine 
carbon sources and assimilation pathways of species in microbial communities. Proceedings of the National Academy 
of Sciences 115 (24), E5576-E5584. 
doi: [https://doi.org/10.1073/pnas.1722325115 ](https://doi.org/10.1073/pnas.1722325115 )

Kleiner M, Kouris A, Jensen M, Liu Y, McCalder J, Strous M (2023) Ultra-sensitive Protein-SIP to quantify activity
and substrate uptake in microbiomes with stable isotopes. Microbiome 11, 24.
doi: [https://doi.org/10.1186/s40168-022-01454-1](https://doi.org/10.1186/s40168-022-01454-1)

M Kösters, J Leufken, S Schulze, K Sugimoto, J Klein, R P Zahedi, M Hippler, S A Leidel, C Fufezan; pymzML v2.0: 
introducing a highly compressed and seekable gzip format, Bioinformatics, 
doi: [https://doi.org/10.1093/bioinformatics/bty046](https://doi.org/10.1093/bioinformatics/bty046)
