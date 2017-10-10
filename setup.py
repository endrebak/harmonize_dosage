from setuptools import setup

from harmonize.version import __version__

setup(
    name = "harmonize",
    packages = ["harmonize"],
    scripts = ["bin/harmonize"],
    version = __version__,
    description = "Harmonize exposure and outcome for Mendelian Randomization.",
    author = "Endre Bakken Stovner",
    author_email = "endrebak85@gmail.com",
    url = "http://github.com/endrebak/harmonize",
    keywords = ["snp", "gwas", "strand"],
    license = ["GPL-3.0"],
    install_requires = ["pandas>=0.16"],
    classifiers = [
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"],
    long_description = ("Harmonize exposure and outcome for Mendelian Randomization. \n"
                        "See the URL for examples and docs.")
)
