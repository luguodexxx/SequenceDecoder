#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/4/10 10:17 PM
__author__ = 'Zhou Ran'

from setuptools import setup, find_packages

setup(name='dlaa',
      version='1.0.0',
      python_requires='>3.6',
      description='DL_Feature',
      author='Ran zhou',
      author_email='ranzhou1005@gmail.com',
      url='https://github.com/luguodexxx/DL_DNA',
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      keywords='DL_feature',
      packages=find_packages(),
      install_requires=[
          'pandas',
          'numpy'
      ],
      entry_points={
          'console_scripts': [
              'dlaa=DLAA.cli:cli'
          ],
      },
      )
