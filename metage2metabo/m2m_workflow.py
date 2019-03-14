#!/usr/bin/env python
# Copyright (c) 2019, Clemence Frioux <clemence.frioux@inria.fr>
#
# This file is part of metage2metabo.
#
# metage2metabo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# metage2metabo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with metage2metabo.  If not, see <http://www.gnu.org/licenses/>.
# -*- coding: utf-8 -*-

import argparse
import time

from metage2metabo import utils
from os import listdir
from os.path import isfile, join


###############################################################################
#
message = """
Description here
"""

pusage = """
Usage here
"""

requires = """
Requirements here"
"""
#
###############################################################################


def run_workflow():
    """description
    """
    parser = argparse.ArgumentParser(description=message, usage=pusage, epilog=requires)
    parser.add_argument("-g", "--genomes",
                        help="annotated genomes repository", required=True)

    args = parser.parse_args()

    print('Hello world, I do nothing more so far ¯\_(ツ)_/¯')
