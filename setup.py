from setuptools import setup, find_packages
setup(
    name = "planck",
    version = "0.3",
    packages = ['planck'], 

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires = ['docutils>=0.3'],

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
    },

    # metadata for upload to PyPI
    author = "Andrea Zonca",
    author_email = "code@andreazonca.com",
    description = "Python package for working with Planck data",
    license = "PSF",
    keywords = "Planck science data",
    url = "http://andreazonca.com/software/planck/",   # project home page, if any

    # could also include long_description, download_url, classifiers, etc.
)
