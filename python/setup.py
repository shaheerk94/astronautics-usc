from setuptools import setup, find_packages

setup(
    name="astronautics-usc",
    version="0.1.0",
    description="Collection of space and orbital mechanics functions for ASTE 520 coursework",
    author="Shaheer Khan",
    packages=find_packages(where=".", exclude=["test"]),
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "matplotlib>=3.4",
        "pint"
    ],
    python_requires=">=3.8",
)
