import sys
import os
from setuptools import setup

main_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(main_dir, "GBpy"))
import GBpy
# del sys.path[0]


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='GBpy',
    version='0.1.0',
    author='Arash Dehghan Banadaki, Srikanth Patala',
    author_email='adehgha@ncsu.edu, spatala@ncsu.edu',
    description="GBpy is an opensource python package for calculating the geometric properties of interfaces in crystals.",
    long_description=read('README.rst'),
    url='https://github.com/adehgha/GBpy',
    download_url = 'https://github.com/adehgha/GBpy/tarball/0.1.0',
    platforms='any',
    requires=['numpy', 'scipy'],
    keywords = ['bicrystallography', 'interfaces', 'grain boundaries'],
    classifiers=['Development Status :: 2 - Pre-Alpha', 'Topic :: Utilities'],
    license='License :: GNU-GPL',
)

# from distutils.core import setup
# with open('README.rst') as f:
#     long_description = f.read()

# setup(name='GBpy',
#     version='0.1.0',
#     author='Arash Dehghan Banadaki, Srikanth Patala',
#     author_email='adehgha@ncsu.edu, spatala@ncsu.edu',
#     description="GBpy is an opensource python package for calculating the geometric properties of interfaces in crystals.",
#     long_description=long_description,
#     url='https://github.com/adehgha/GBpy',
#     download_url = 'https://github.com/adehgha/GBpy/tarball/0.1.0',
#     platforms='any',
#     requires=['numpy', 'scipy'],
#     keywords = ['bicrystallography', 'interfaces', 'grain boundaries'],
#     classifiers=['Development Status :: 2 - Pre-Alpha', 'Topic :: Utilities'],
#     license='License :: GNU-GPL',
# )
