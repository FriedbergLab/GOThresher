# flake8: noqa
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
# CreateTime: 2017-06-05 15:13:04

#! /usr/bin/env python
setup(name='crytto',
      description='Remove bias from GAF files based on annotation information content, GO evidence, annotation source, number of proteins annotated from a given source, an date',
      long_description=long_description,
      version='0.1',
      url='https://github.com/Rinoahu/debias',
      author='idoerg',
      author_email='idoerg@gmail.com',
      license='GPL',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: System Administrators',
          'License :: OSI Approved :: GPL License',
          'Programming Language :: Python :: 2.7'

      ],
      packages=['crytto'],
      install_requires=[
          'PyYAML>=3.11',
          'sh>=1.11'
      ],
      entry_points={
          'console_scripts': [
              'encrypt=crytto.main:run'
          ]
      }
)
