#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

setup(
    name='QEplayground',
    version='0.1.0',
    description="Quantum Espresso playground",
    long_description=readme + '\n\n' + history,
    author="Claudio Attaccalite and Elena Cannuccia",
    author_email='claudio.attaccalite@gmail.com',
    url='https://github.com/attacc/QEplayground',
    packages=find_packages(include=['QEplayground']),
    entry_points={
        'console_scripts': [
            'QEplayground=QEplayground.cli:main'
        ]
    },
    include_package_data=True,
    install_requires=['numpy'],
    license="GNU General Public License v3",
    zip_safe=False,
    keywords='QEplayground',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    python_requires='>=3',
)
