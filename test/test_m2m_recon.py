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


def tca_reactions():
    return ["RXN-14971", "MALATE-DEH-RXN", "ISOCITDEH-RXN", "MALATE-DEHYDROGENASE-ACCEPTOR-RXN",
                    "ACONITATEDEHYDR-RXN", "CITSYN-RXN", "ACONITATEHYDR-RXN", "2OXOGLUTARATEDEH-RXN", "SUCCCOASYN-RXN", "FUMHYDR-RXN"]


def fabo_reactions():
    return ["ACYLCOASYN-RXN", "ACYLCOADEHYDROG-RXN", "ENOYL-COA-DELTA-ISOM-RXN", "ENOYL-COA-HYDRAT-RXN",
                    "OHBUTYRYL-COA-EPIM-RXN", "OHACYL-COA-DEHYDROG-RXN", "KETOACYLCOATHIOL-RXN"]


def test_m2m_recon_call():
    """
    Test m2m recon when called in terminal.
    """
    subprocess.call(['mpwt', '--delete', 'fatty_acid_beta_oxydation_icyc,tca_cycle_ecolicyc'])
    subprocess.call(['m2m', 'recon', '-g', 'recon_data', '-o', 'recon_data_output', '-c', '1'])

    reader = SBMLReader()
    document = reader.readSBML('recon_data_output/sbml/tca_cycle_ecoli.sbml')
    expected_tca_reactions = [convert_from_coded_id(reaction.getId())[0] for reaction in document.getModel().getListOfReactions()]
    assert set(tca_reactions()).issubset(set(expected_tca_reactions))

    reader = SBMLReader()
    document = reader.readSBML('recon_data_output/sbml/fatty_acid_beta_oxydation_I.sbml')
    expected_fabo_reactions = [convert_from_coded_id(reaction.getId())[0] for reaction in document.getModel().getListOfReactions()]
    assert set(fabo_reactions()).issubset(set(expected_fabo_reactions))

    shutil.rmtree('recon_data_output')
    subprocess.call(['mpwt', '--delete', 'fatty_acid_beta_oxydation_icyc,tca_cycle_ecolicyc'])

