#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2019-2023 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

from setuptools import setup

setup(
    name='Metage2Metabo',
    url='https://github.com/aureme/metage2metabo',
    license='GPLv3+',
    description=
    'Automatic reconstruction of draft metabolic networks with Pathway Tools and graph-based metabolic analysis',
    long_description=
    'metage2metabo is a Python3 workflow to perform graph-based metabolic analysis starting from annotated genomes. It uses Pathway Tools in a automatic and parallel way to reconstruct metabolic networks for a large number of genomes. The obtained metabolic networks are then analyzed individually and collectively in order to get the added value of cooperation in microbiota over individual metabolism, and to identify and screen interesting organisms among all. \
More information on usage and troubleshooting on Github: https://github.com/aureme/metage2metabo',
    author='AuReMe',
    author_email='gem-aureme@inria.fr',
    packages=['metage2metabo', 'metage2metabo.m2m', 'metage2metabo.m2m_analysis'],
    package_dir={'metage2metabo': 'metage2metabo', 'metage2metabo.m2m': 'metage2metabo/m2m' , 'metage2metabo.m2m_analysis':'metage2metabo/m2m_analysis'},
    package_data={'metage2metabo': ['workflow_data/workflow_genomes.tar.gz',
                                    'workflow_data/seeds_workflow.sbml']},
    entry_points={
        'console_scripts': [
            'm2m = metage2metabo.__main__:main',
            'm2m_analysis = metage2metabo.__main_analysis__:main',
        ]
    },
    install_requires=['miscoto', 'menetools', 'mpwt', 'padmet'],
)
