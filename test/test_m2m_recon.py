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

def fabo_reactions():
    return ["ACYLCOASYN-RXN", "ACYLCOADEHYDROG-RXN", "ENOYL-COA-DELTA-ISOM-RXN", "ENOYL-COA-HYDRAT-RXN",
                    "OHBUTYRYL-COA-EPIM-RXN", "OHACYL-COA-DEHYDROG-RXN", "KETOACYLCOATHIOL-RXN"]


def test_m2m_recon_call():
    """
    Test m2m recon when called in terminal.
    """
    subprocess.call(['mpwt', '--delete', 'fatty_acid_beta_oxydation_icyc'])
    subprocess.call(['m2m', 'recon', '-g', 'recon_data', '-o', 'recon_data_output', '-c', '1', '-p'])

    reader = SBMLReader()
    document = reader.readSBML('recon_data_output/sbml/fatty_acid_beta_oxydation_I.sbml')
    expected_fabo_reactions = [convert_from_coded_id(reaction.getId())[0] for reaction in document.getModel().getListOfReactions()]
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    padmet = PadmetSpec('recon_data_output/padmet/fatty_acid_beta_oxydation_I.padmet')
    fabo_rxns = [node.id for node in padmet.dicOfNode.values() if node.type == "reaction"]
    assert set(fabo_reactions()).issubset(set(fabo_rxns))

    shutil.rmtree('recon_data_output')
    subprocess.call(['mpwt', '--delete', 'fatty_acid_beta_oxydation_icyc'])

