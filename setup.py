import sys
import os
from setuptools import setup

main_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(main_dir, "GBpy"))
import GBpy
# del sys.path[0]

setup(name='GBpy',
    version='1.0',
    description="GBpy is an opensource python package for calculating the geometric properties of interfaces in crystals.",
    author='Arash Dehghan Banadaki, Srikanth Patala',
    author_email='adehgha@ncsu.edu, spatala@ncsu.edu',
    license='Copyright the authors.'
)
