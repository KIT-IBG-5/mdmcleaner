import pathlib
import sys
from setuptools import find_packages, setup
from mdmcleaner import _version

cwd = pathlib.Path(__file__).parent
readme = (cwd / "README.md").read_text()

setup(
    name = "mdmcleaner",
    version = _version,
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
    #package_dir = {"mdmcleaner" : "lib"},
    packages = ["mdmcleaner"],
    #scripts = ["mdmcleaner"],
    include_package_data=True,
    license = "GPL-3.0",
    entry_points={"console_scripts":["mdmcleaner=mdmcleaner.mdmcleaner:main"]},
)
