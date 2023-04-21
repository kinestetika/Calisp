import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='calisp',
    version='3.0.11',
    packages=setuptools.find_packages(where='src'),
    url='https://github.com/kinestetika/Calisp',
    license='MIT',
    author='Marc Strous',
    author_email='mstrous@ucalgary.ca',
    description='Isotope analysis of proteomics data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Natural Language :: English',
                 'Operating System :: OS Independent',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.9',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='proteomics isotope mass-spectrometry 13C',
    project_urls={'Source': 'https://github.com/kinestetika/Calisp'},
    package_dir={'': 'src'},
    python_requires='>=3.6',
    install_requires=['numpy', 'scipy', 'pandas', 'tqdm', 'pymzml', 'pyarrow'],
    license_files=('LICENSE.txt',),
    extras_require={  # Optional
        'dev': ['setuptools', 'build', 'twine'],
        'test': ['jupyter', 'matplotlib', 'jinja2'],
    },
    entry_points={  # Optional
        'console_scripts': [
            'calisp=calisp.main:main',
        ],
    }
)
