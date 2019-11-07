import sys
import os
from setuptools import setup, find_packages

main_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(main_dir, "GBpy"))

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

setup(name='GBpy',
    version='0.1.2',
    author='Arash Dehghan Banadaki, Srikanth Patala',
    author_email='adehgha@ncsu.edu, spatala@ncsu.edu',
    description="GBpy is an opensource python package for calculating the geometric properties of interfaces in crystals.",
    long_description=read('README.rst'),
    url='https://github.com/adehgha/GBpy',
    download_url = 'https://github.com/adehgha/GBpy/tarball/0.1.2',
    platforms='any',
    requires=['numpy', 'scipy'],
    packages = find_packages(),
    keywords = ['bicrystallography', 'interfaces', 'grain boundaries'],
    classifiers=['Development Status :: 2 - Pre-Alpha', 'Topic :: Utilities'],
    license='License :: GNU-GPL',
    cmdclass = {'build_py' : build_py}
)
