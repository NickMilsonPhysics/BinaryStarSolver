from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()


setup(
    name="BinaryStarSolver",
    version="0.0.9",
    author="Nicholas Milson",
    author_email="nick.milson@dal.ca",
    description="Solves for the orbital elements of binary stars, given radial velocity time series",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/NickMilsonPhysics/binarystarsolve",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
