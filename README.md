## Calisp.py, version 3.0

Calisp.py (Calgary approach to isotopes in proteomics) is a program that estimates isotopic composition (e.g. 13C/12C,
delta13C, 15N/14N etc) of peptides from proteomics mass spectrometry data. Input data consist of mzML files and 
files with peptide spectrum matches.

Calisp was originally developed in Java. This newer, python version is much more concise, it consists of only six .py 
files, ~1,000 lines of code, compared to about fifty files and ~10,000 lines of code for the Java program. In addition,
the Java program depends on another program, mcl, whereas calisp.py is a purely python program, which is easy to
install, see below. The conciseness of the code makes the python version more transparent, easier to maintain and 
easier to further develop. I consider Calisp.py is the successor of the Java version. One of the reasons to shift to 
python is the possibility to more effectively develop machine learning approaches to filter out noisy spectra. Future
work will explore that possibility.

Benchmarking of Calisp.py has been completed. It works well. Parsing of .mzid files still needs to be implemented. 

Calisp.py depends on numpy, scipy, pandas, tqdm, pymzml, pyarrow. It outputs the data as a Pandas DataFrame saved in 
(binary) feather format. From there, the user can explore and visualize the results in a Jupyter notebook. For that, a
tutorial wil be provided.

Compared to previous versions of calisp, the workflow has been simplified. Calisp.py does not filter out any spectra, or
adds up spectra to reduce noise - like the Java version does. It simply estimates the ratio for the target isotopes 
(e.g. 13C/12C) for every spectrum it can find. It estimates this ratio based on neutron abundance and using fast fourier
transforms. The former applies to stable isotope probing experiments. The latter applies to natural abundances, or to 
isotope probing experiments with very little added label (e.g. using substrates wiht only 1% additional 13C). The
motivation for omitting filtering is that collecting all spectra, including bad ones, will enable training of machine 
learning classifiers. Also, because it was shown that the median provides better estimates for species in microbial
communities than the mean, adding up spectra to improve precision has lost its purpose. There is more power 
(and sensitivity) in numbers.

Because no data are filtered out and no spectra get added up, calisp.py, analyzes at least ten times as many spectra
compared to the Java version. That means calisp.py is about ten times slower, it takes about 5-10 min per .mzML file on 
a Desktop computer. The user can perform filtering of the data using the Jupyter notebook as desired. For natural 
abundance data, it works well to only use those spectra that have a FFT fitting error ("error_fft") of less than 
0.001, as previously shown for the java program.

**Installation:**

>python -m virtualenv calisp

>source calisp/bin/activate

>pip install calisp

**Usage:**

>source calisp/bin/activate

>calisp.py --spectrumFile [path to .mzML file or folder with .mzML files] --peptideFile [path to .mzML file or folder 
 with .mzML files] --outputFile [path to where calisp.py will save results files] --threads [# of threads used, 
 default 4] --isotope [15N, 3H etc, default 13C], --organismDelimiter [character that separates an organism ID from the
 remainder of protein IDs, default '_'], --massAccuracy [accuracy of peak m/z identifications, default 10 ppm]
 --compute_clumps [use if you want to compute clumpiness]

**Column names of the Pandas DataFrame created by calisp.py:**

In the saved dataframe, each row contains one spectrum, defined by the following columns:

 'experiment', 'ms_run', 'bins', 'proteins', 'peptide', 'peptide_mass', 'C', 'N', 'O', 'H', 'S',
 'psm_id','psm_mz', 'psm_charge', 'psm_neutrons', 'psm_rank', 'psm_precursor_id',
 'psm_precursor_mz', 'spectrum_charge', 'spectrum_precursor_id', 'spectrum_total_intensity',
 'spectrum_peak_count', 'spectrum_median_peak_spacing', 'spectrum_mass_irregularity',
 'ratio_na', 'ratio_fft','error_fft', 'error_clumpy'
 'flag_peptide_contains_sulfur', 'flag_peptide_has_modifications',
 'flag_peptide_assigned_to_multiple_bins', 'flag_peptide_assigned_to_multiple_proteins',
 'flag_peptide_mass_and_elements_undefined', 'flag_psm_has_low_confidence','flag_psm_is_ambiguous',
 'flag_spectrum_is_contaminated', 'flag_spectrum_is_wobbly', 'flag_peak_at_minus_one_pos'
 'i0', ..., 'i19', 'm0', ..., 'm19', 'c1' ... 'c6'

The final 46 columns contain the normalized peak intensities ('i0' ...) and m/z of the spectrum's peaks ('m0' ...), 
as well as the inferred clumpiness of the isotopes ('c1' ...).

**Please cite:**

Kleiner M, Dong X, Hinzke T, Wippler J, Thorson E, Mayer B, Strous M (2018) A metaproteomics method to determine 
carbon sources and assimilation pathways of species in microbial communities. Proceedings of the National Academy 
of Sciences 115 (24), E5576-E5584.
