# -*- coding: utf-8 -*-

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# Get the long description from the README file
import os
import codecs
here = os.path.abspath(os.path.dirname(__file__))
with codecs.open('README.rst', encoding='utf-8') as f:
	long_description = f.read()



# setup arguments
setup(
	name='fixif.FXPF',
	version='1.0',
	description='Python Wrapper for the Fixed-Point Formats (FXPF) library',
	long_description=long_description,
	url='https://github.com/fixif/FXPF',
	author='A. Volkova',
	author_email='anastasia.volkova@inria.fr',
	classifiers=[
		'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Topic :: Scientific/Engineering',
		'Topic :: Software Development :: Embedded Systems',
		'Topic :: Software Development :: Code Generators',
		'Programming Language :: Python :: 2',
	],
	keywords='fixed-point formats, FXPF',
	packages=find_packages(exclude=['tests']),
	install_requires=['pytest', 'pytest-cov', 'coveralls', 'numpy', 'scipy'],   # packages only used for testing
	project_urls={
		'Bug Reports': 'https://github.com/FiXiF/FXPF/issues',
		'Source': 'https://github.com/FiXiF/FXPF/',
	},
)
