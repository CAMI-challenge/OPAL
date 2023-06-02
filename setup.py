from setuptools import setup, find_packages
from version import __version__
import glob


def dependencies():
    with open('requirements.txt', 'r') as f:
        return f.read().splitlines()

setup(
    name             = 'cami-opal',
    version          = __version__,
    description      = 'OPAL: Open-community Profiling Assessment tooL',
    author           = 'CAMI',
    author_email     = 'support@cami-challenge.org',
    url              = 'http://cami-challenge.org',
    license          = 'Apache-2.0 License',
    scripts          = glob.glob('*.py'),
    install_requires = dependencies(),
    packages         = find_packages(),
    classifiers = [
        'Natural Language :: English',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX'
    ]
)
