#!/usr/bin/env python3
"""
Palace: 3D Finite Element Solver for Computational Electromagnetics

Setup script for building and installing Palace as a Python package.
"""

from setuptools import setup, find_packages
import os
import re

def read(fname):
    """Read file contents."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def get_version():
    """Extract version from CMakeLists.txt."""
    with open('CMakeLists.txt', 'r') as f:
        content = f.read()

    # Look for version pattern in CMakeLists.txt
    version_match = re.search(r'project\s*\(\s*[Pp]alace\s+VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)', content)
    if version_match:
        return version_match.group(1)

    # Fallback version if not found
    return "0.1.0"

# Read the long description from README
long_description = read('README.md')

setup(
    name="palace-fem",
    version=get_version(),
    author="AWS Center for Quantum Computing",
    author_email="palace-maint@amazon.com",
    description="3D Finite Element Solver for Computational Electromagnetics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/awslabs/palace",
    project_urls={
        "Bug Tracker": "https://github.com/awslabs/palace/issues",
        "Documentation": "https://awslabs.github.io/palace/",
        "Source Code": "https://github.com/awslabs/palace",
    },
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme",
        ],
        "visualization": [
            "matplotlib>=3.3.0",
            "paraview>=5.9.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "palace=palace.cli:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="finite element electromagnetics simulation physics",
)
