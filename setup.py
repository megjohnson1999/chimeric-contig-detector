from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="chimeric-detective",
    version="1.0.0",
    author="Chimeric Detective Team",
    description="A comprehensive tool for detecting and resolving chimeric contigs in viral metagenomic assemblies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/megjohnson1999/chimeric-contig-detector",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79",
        "pysam>=0.19.0",
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "plotly>=5.0.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "click>=8.0.0",
        "tqdm>=4.62.0",
        "jinja2>=3.0.0",
    ],
    entry_points={
        "console_scripts": [
            "chimeric_detective=chimeric_detective.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "chimeric_detective": ["templates/*.html", "templates/*.css", "templates/*.js"],
    },
)