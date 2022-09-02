import sys
from setuptools import setup

if sys.version_info < (3, 8):
    print("WARNING: Python 3.8 or newer is highly recommended. If you encounter errors, try updating python")

install_requires = [
    "cython>=0.29.30, <3.0.0",
    "numpy>=1.20.1",
    "matplotlib>=3.3.4",
    "scipy==1.6.1",
    "liftover==1.1.13",
    "python-magic==0.4.24",
    "requests>=2.22.0",
    "tqdm>=4.60.0, <5.0.0",
]

classifiers = """
Programming Language :: Python :: 3.8
Programming Language :: Bash :: 4
Programming Language :: Awk
Development Status :: 5 - Production/Stable 
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
"""

if __name__ == "__main__":
    setup(
        name="SumStatsRehab",
        version="1.2.1",
        description="GWAS summary statistics files QC tool",
        url="https://github.com/Kukuster/SumStatsRehab",
        license="MIT",
        classifiers=classifiers.split("\n"),
        zip_safe=False,
        py_modules=['SumStatsRehab',
            'lib/env',
            'lib/file',
            'lib/loop_fix',
            'lib/math_utils',
            'lib/prepare_GWASSS_columns',
            'lib/prepare_two_dbSNPs',
            'lib/report_utils',
            'lib/sort_GWASSS_by_ChrBP',
            'lib/sort_GWASSS_by_rsID',
            'lib/standard_column_order',
            'lib/utils',
            'lib/validate_GWASSS_entries',
        ],
        install_requires=install_requires,
        entry_points={"console_scripts": ["SumStatsRehab=SumStatsRehab:main"]},
    )
