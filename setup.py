import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='calisp',
    version='0.0',
    packages=setuptools.find_packages(where='src'),
    url='https://github.com/kinestetika/Calisp',
    license='MIT',
    author='Marc Strous',
    author_email='mstrous@ucalgary.ca',
    description='Isotope analysis of proteomics data',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    classifiers=['Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'Framework :: Jupyter',
                 'Natural Language :: English',
                 'Operating System :: OS Independent',
                 'Intended Audience :: Scientists',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.9',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='proteomics isotope mass-spectrometry 13C',
    project_urls={'Source': 'https://github.com/kinestetika/Calisp'},
    package_dir={'': 'src'},
    python_requires='>=3.6',
    install_requires=['numpy', 'scipy', 'pandas', 'tqdm', 'pymzml', 'pyarrow'],
)