import pathlib
import sys
from setuptools import find_packages, setup
import mdmcleaner

cwd = pathlib.Path(__file__).parent
readme = (cwd / "README.md").read_text()

setup(
    name = "mdmcleaner",
    version = mdmcleaner.__version__,
    description = "A pipeline for the assessment, classification and refinement of microbial dark matter SAGs and MAGs",
	long_description = readme,
	long_description_content_type="text/markdown",
	url = "https://github.com/KIT-IBG-5/mdmcleaner",
    author = "John Vollmers",
    author_email = "john.vollmers@kit.edu",
    install_requires=["biopython", "numpy"],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    package_dir = {"mdmcleaner" : "mdmcleaner"},
    packages = ["mdmcleaner", "mdmcleaner/hmms", "mdmcleaner/hmms/arch", "mdmcleaner/hmms/prok", "mdmcleaner/hmms/bact"],
    package_data = {"mdmcleaner" : ["hmms/cutofftable_combined.tsv", "hmms/arch/*.hmm", "hmms/prok/*.hmm", "hmms/bact/*.hmm"]},
    #scripts = ["mdmcleaner"],
    include_package_data=True,
    license = "GPL-3.0",
    entry_points={"console_scripts":["mdmcleaner=mdmcleaner.mdmcleaner:main"]},
)
