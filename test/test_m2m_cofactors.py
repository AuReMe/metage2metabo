#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test m2m recon on genbank files containing E. coli genes implied in the TCA cycle and in the Fatty Acid Beta oxydation.
Need an environment with Pathway-Tools installed.
"""

import os
import shutil
import subprocess

from libsbml import SBMLReader
from padmet.utils.sbmlPlugin import convert_from_coded_id
from padmet.classes import PadmetSpec


def test_m2m_cofactors_call():
    """
    Test m2m recon when called in terminal.
    """
    expected_reactions = ['reaction__cof__', 'R_reaction']
    expedcted_metabolites = ['M_ADP_c', 'M_ATP_c', 'M_metabolite_A_c', 'M_metabolite_B_c', 'M_PROTON_c', 'ATP_c__cof__', 'ADP_c__cof__']

    input_folder = os.path.join('metabolic_data', 'sbml_cofactors')
    cofactor_file = os.path.join('metabolic_data', 'cofactors.tsv')

    subprocess.call(['m2m', 'cofactors', '-n', input_folder, '--cofactors', cofactor_file, '-o', 'cofactors_output', '-c', '1'])

    sbml_output_folder = os.path.join('cofactors_output', 'mocked_sbmls')
    sbml_output_file = os.path.join(sbml_output_folder, 'test.sbml')

    reader = SBMLReader()
    document = reader.readSBML(sbml_output_file)
    found_reactions = [reaction.getId() for reaction in document.getModel().getListOfReactions()]

    assert set(expected_reactions) == set(found_reactions)
    found_metabolites = [reaction.getId() for reaction in document.getModel().getListOfSpecies()]
    assert set(expedcted_metabolites) == set(found_metabolites)

    shutil.rmtree('cofactors_output')
