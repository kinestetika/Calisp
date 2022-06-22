## Calisp.py, version 3.0.8

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
 --compute_clumps [use only if you want to compute clumpiness]

**Column names of the Pandas DataFrame created by calisp.py:**

In the saved dataframe, each row contains one isotopic pattern, defined by the following columns:

 'experiment', 'ms_run', 'bins', 'proteins', 'peptide', 'peptide_mass', 'C', 'N', 'O', 'H', 'S',
 'psm_id','psm_mz', 'psm_charge', 'psm_neutrons', 'psm_rank', 'psm_precursor_id',
 'psm_precursor_mz', 'spectrum_charge', 'spectrum_precursor_id', 'spectrum_total_intensity',
 'spectrum_peak_count', 'spectrum_median_peak_spacing', 'spectrum_mass_irregularity',
 'ratio_na', 'ratio_fft','error_fft', 'error_clumpy'
 'flag_peptide_contains_sulfur', 'flag_peptide_has_modifications',
 'flag_peptide_assigned_to_multiple_bins', 'flag_peptide_assigned_to_multiple_proteins',
 'flag_peptide_mass_and_elements_undefined', 'flag_psm_has_low_confidence','flag_psm_is_ambiguous',
 'flag_pattern_is_contaminated', 'flag_pattern_is_wobbly', 'flag_peak_at_minus_one_pos'
 'i0', ..., 'i19', 'm0', ..., 'm19', 'c1' ... 'c6'

The final 46 columns contain the normalized peak intensities ('i0' ...) and m/z of the patterns' peaks ('m0' ...), 
as well as the inferred clumpiness of the isotopes ('c1' ...).

calisp.py was developed using [PyCharm comunity edition](https://www.jetbrains.com/pycharm/).

**Please cite:**

Kleiner M, Dong X, Hinzke T, Wippler J, Thorson E, Mayer B, Strous M (2018) A metaproteomics method to determine 
carbon sources and assimilation pathways of species in microbial communities. Proceedings of the National Academy 
of Sciences 115 (24), E5576-E5584. 
doi: [https://doi.org/10.1073/pnas.1722325115 ](https://doi.org/10.1073/pnas.1722325115 )

Kleiner M, Kouris A, Jensen M, Liu Y, McCalder J, Strous M (2021) Ultra-sensitive Protein-SIP to quantify activity 
and substrate uptake in microbiomes with stable isotopes. bioRxiv.
doi: [https://doi.org/10.1101/2021.03.29.437612](https://doi.org/10.1101/2021.03.29.437612)

M Kösters, J Leufken, S Schulze, K Sugimoto, J Klein, R P Zahedi, M Hippler, S A Leidel, C Fufezan; pymzML v2.0: 
introducing a highly compressed and seekable gzip format, Bioinformatics, 
doi: [https://doi.org/10.1093/bioinformatics/bty046](https://doi.org/10.1093/bioinformatics/bty046)
