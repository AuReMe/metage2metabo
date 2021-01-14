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

import logging
import os
import time

from metage2metabo.m2m_analysis.enumeration import enumeration_analysis
from metage2metabo.m2m_analysis.solution_graph import graph_analysis
from metage2metabo.m2m_analysis.graph_compression import powergraph_analysis

logger = logging.getLogger(__name__)


def run_analysis_workflow(sbml_folder, target_folder_file, seed_file, output_dir, taxon_file, oog_jar, host_file=None, taxonomy_level="phylum"):
    """Run the whole m2m_analysis workflow

    Args:
        sbml_folder (str): sbml directory
        target_folder_file (str): targets file or folder containing multiple sbmls
        seed_file (str): seeds file
        output_dir (str): results directory
        taxon_file (str): mpwt taxon file for species in sbml folder
        oog_jar (str): path to OOG jar file
        host_file (str): metabolic network file for host
    """
    starttime = time.time()

    json_file_folder = enumeration_analysis(sbml_folder, target_folder_file, seed_file, output_dir, host_file)

    gml_output = graph_analysis(json_file_folder, target_folder_file, output_dir, taxon_file, taxonomy_level)

    powergraph_analysis(gml_output, output_dir, oog_jar, taxon_file, taxonomy_level)

    logger.info(
        "--- m2m_analysis runtime %.2f seconds ---\n" % (time.time() - starttime))
