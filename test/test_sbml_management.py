#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test sbml management functions used by m2m.
"""

import os
import pytest
from libsbml import SBMLReader
from metage2metabo.sbml_management import create_species_sbml


def test_create_species_sbml():
    metabolites = set(['M_A_c', 'M_B_c', 'M_C_c', 'M_D_c'])
    metabolite_file = os.path.join('test_metabolite.sbml')

    create_species_sbml(metabolites, metabolite_file)

    reader = SBMLReader()
    document = reader.readSBML(metabolite_file)
    created_metabolites = set([specie.getId() for specie in document.getModel().getListOfSpecies()])

    assert created_metabolites == metabolites

    os.remove(metabolite_file)


def test_create_species_sbml_invalid_id():
    metabolites = set(['M_A+_c', 'M_B_\c', 'M_C!=_c', 'M_D_c'])
    metabolite_file = os.path.join('test_metabolite.sbml')

    with pytest.raises(SystemExit) as pytest_exit:
        create_species_sbml(metabolites, metabolite_file)

    assert pytest_exit.type == SystemExit
    assert pytest_exit.value.code == 1


if __name__ == "__main__":
    test_create_species_sbml()
