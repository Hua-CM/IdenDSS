# -*- coding: utf-8 -*-
# @Time : 2022/9/25 21:00
# @Author : Zhongyi Hua
# @FileName: setup.py.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

__version__ = "0.2.4"

setuptools.setup(
    name="IdenDSS",
    version=__version__,
    author="Zhongyi Hua",
    author_email="njbxhzy@hotmail.com",
    description="DNA Signature Sequence (DSS) Identification and Application Software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Hua-CM/IdenDSS",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "biopython>=1.78",
        "pandas>=1.0.0",
        "numpy>=1.20.0"],
    include_package_data=True,
    entry_points={'console_scripts': ['IdenDSS = IdenDSS.IDSS:main']},
    package_dir={'': "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.8",
)
