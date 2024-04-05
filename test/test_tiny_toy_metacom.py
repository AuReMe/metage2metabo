#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test m2m metacom on a tiny dataset that is easily visualized and understood.
"""

import os
import shutil
import json
import sys
from libsbml import SBMLReader
from metage2metabo.m2m.m2m_workflow import metacom_analysis

ISCOPE = {
    "bact1": [
        "M_A_c",
        "M_B_c",
        "M_R_c",
        "M_H_c",
        "M_D_c"
    ],
    "bact6": [],
    "bact4": [
        "M_B_c"
    ],
    "bact5": [],
    "bact2": [
        "M_B_c",
        "M_N_c",
        "M_A_c",
        "M_S2_c"
    ],
    "bact3": [
        "M_A_c",
        "M_B_c"
    ]
}

PRODUCED_SEEDS_ISCOPE = {
    "individually_producible_seeds": {
        "bact1": [],
        "bact2": [
            "M_S2_c"
        ],
        "bact3": [],
        "bact4": [],
        "bact5": [],
        "bact6": []
    },
    "seeds_absent_in_metabolic_network": {
        "bact1": [],
        "bact2": [],
        "bact3": [],
        "bact4": [
            "M_S1_c"
        ],
        "bact5": [
            "M_S2_c"
        ],
        "bact6": [
            "M_S2_c"
        ]
    }
}

REV_ISCOPE = {
    "M_A_c": [
        "bact1",
        "bact2",
        "bact3"
    ],
    "M_B_c": [
        "bact1",
        "bact4",
        "bact2",
        "bact3"
    ],
    "M_D_c": [
        "bact1"
    ],
    "M_H_c": [
        "bact1"
    ],
    "M_N_c": [
        "bact2"
    ],
    "M_R_c": [
        "bact1"
    ],
    "M_S2_c": [
        "bact2"
    ]
}

COMSCOPE = {
    "com_scope": [
        "M_X_c",
        "M_N_c",
        "M_A_c",
        "M_K_c",
        "M_F_c",
        "M_B_c",
        "M_C_c",
        "M_R_c",
        "M_H_c",
        "M_V_c",
        "M_S1_c",
        "M_D_c",
        "M_E_c",
        "M_G_c",
        "M_S2_c"
    ]
}

REV_CSCOPE = {
    "M_A_c": [
        "bact1",
        "bact2",
        "bact3"
        ],
    "M_B_c": [
        "bact1",
        "bact2",
        "bact3",
        "bact4"
        ],
    "M_C_c": [
        "bact5",
        "bact6"
        ],
    "M_D_c": [
        "bact1"
        ],
    "M_E_c": [
        "bact2",
        "bact3"
        ],
    "M_F_c": [
        "bact1"
        ],
    "M_G_c": [
        "bact2",
        "bact3"
        ],
    "M_H_c": [
        "bact1",
        "bact2",
        "bact3"
        ],
    "M_K_c": [
        "bact5",
        "bact6"
        ],
    "M_N_c": [
        "bact2",
        "bact3"
        ],
    "M_R_c": [
        "bact1"
        ],
    "M_S1_c": [
        "bact3"
        ],
    "M_S2_c": [
        "bact2"
        ],
    "M_V_c": [
        "bact5"
        ],
    "M_X_c": [
        "bact4",
        "bact6"
        ]
}

MICROBE_CONTRIBUTIONS = {
    "bact1": {
        "community_metabolic_gain": [
            "M_F_c"
        ],
        "produced_alone": [
            "M_A_c",
            "M_H_c",
            "M_D_c",
            "M_B_c",
            "M_R_c"
        ],
        "produced_in_community": [
            "M_H_c",
            "M_D_c",
            "M_R_c",
            "M_A_c",
            "M_F_c",
            "M_B_c"
        ]
    },
    "bact2": {
        "community_metabolic_gain": [
            "M_H_c",
            "M_E_c",
            "M_G_c"
        ],
        "produced_alone": [
            "M_B_c",
            "M_N_c",
            "M_S2_c",
            "M_A_c"
        ],
        "produced_in_community": [
            "M_H_c",
            "M_E_c",
            "M_S2_c",
            "M_A_c",
            "M_B_c",
            "M_N_c",
            "M_G_c"
        ]
    },
    "bact3": {
        "community_metabolic_gain": [
            "M_N_c",
            "M_E_c",
            "M_H_c",
            "M_S1_c",
            "M_G_c"
        ],
        "produced_alone": [
            "M_B_c",
            "M_A_c"
        ],
        "produced_in_community": [
            "M_H_c",
            "M_E_c",
            "M_A_c",
            "M_N_c",
            "M_S1_c",
            "M_B_c",
            "M_G_c"
        ]
    },
    "bact4": {
        "community_metabolic_gain": [
            "M_X_c"
        ],
        "produced_alone": [
            "M_B_c"
        ],
        "produced_in_community": [
            "M_B_c",
            "M_X_c"
        ]
    },
    "bact5": {
        "community_metabolic_gain": [
            "M_K_c",
            "M_V_c",
            "M_C_c"
        ],
        "produced_alone": [],
        "produced_in_community": [
            "M_K_c",
            "M_V_c",
            "M_C_c"
        ]
    },
    "bact6": {
        "community_metabolic_gain": [
            "M_X_c",
            "M_K_c",
            "M_C_c"
        ],
        "produced_alone": [],
        "produced_in_community": [
            "M_X_c",
            "M_K_c",
            "M_C_c"
        ]
    }
}

EXPECTED_TARGETS_ADVAL = {
    'M_X_c',
    'M_K_c',
    'M_F_c',
    'M_C_c',
    'M_E_c',
    'M_V_c',
    'M_S1_c',
    'M_G_c'
}

UNION_MINCOM = {
    'bact1', 'bact2', 'bact3', 
    'bact5', 'bact6'
}
INTERSECTION_MINCOM = {
    'bact1'
}
MIN_SIZE_COM = 3
PROD_TARGETS = {
    'M_C_c',
    'M_F_c',
    'M_H_c'
}

UNPROD_TARGETS = {
    'M_foo_c'
}

PRODUCIBILITY_TARGETS = {
    "unproducible": [
        "M_foo_c"
    ],
    "producible": [
        "M_F_c",
        "M_C_c",
        "M_H_c"
    ],
    "indiv_producible": [
        "M_H_c"
    ],
    "individual_producers": {
        "M_H_c": [
            "bact1"
        ]
    },
    "com_only_producers": {
        "M_C_c": [
            "bact5",
            "bact6"
        ],
        "M_H_c": [
            "bact2",
            "bact3"
        ],
        "M_F_c": [
            "bact1"
        ]
    },
    "mincom_producible": [
        "M_F_c",
        "M_C_c",
        "M_H_c"
    ],
    "key_species": [
        "bact2",
        "bact1",
        "bact5",
        "bact3",
        "bact6"
    ],
    "mincom_union_producers": {
        "M_C_c": [
            "bact5",
            "bact6"
        ],
        "M_H_c": [
            "bact1",
            "bact3",
            "bact2"
        ],
        "M_F_c": [
            "bact1"
        ]
    },
    "mincom_inter_producers": {
        "M_H_c": [
            "bact1"
        ],
        "M_F_c": [
            "bact1"
        ]
    }
}

NUMBER_BACT = 6
SIZE_UNION = 7
SIZE_INTERSECTION = 0
SIZE_CSCOPE = 15


def test_m2m_metacom_tiny_toy():
    """
    Test m2m metacom when called from the API using the tiny toy dataset.
    """
    # RUN THE COMMAND
    inppath = 'metabolic_data/tiny_toy'
    respath = 'tiny_metacom_output'
    networks_path = os.path.join(inppath, 'networks')
    seeds_path = os.path.join(inppath, 'seeds_community.sbml')
    targets_path = os.path.join(inppath, 'targets_community.sbml')

    metacom_analysis(sbml_dir = networks_path, out_dir = respath, seeds = seeds_path, host_mn = None, targets_file = targets_path, cpu_number=1)

    # result files
    iscope_file = os.path.join(*[respath, 'indiv_scopes', 'indiv_scopes.json'])
    rev_iscope_file = os.path.join(*[respath, 'indiv_scopes', 'rev_iscope.json'])
    seeds_iscope_file = os.path.join(*[respath, 'indiv_scopes', 'seeds_in_indiv_scopes.json'])
    cscope_file = os.path.join(*[respath, 'community_analysis', 'comm_scopes.json'])
    rev_cscope_file = os.path.join(*[respath, 'community_analysis', 'rev_cscope.json'])
    contrib_microbes_file = os.path.join(*[respath, 'community_analysis', 'contributions_of_microbes.json'])
    addedvalue_file = os.path.join(*[respath, 'community_analysis', 'addedvalue.json'])
    mincom_file = os.path.join(*[respath, 'community_analysis', 'mincom.json'])
    targets_file = os.path.join(*[respath, 'community_analysis', 'targets.sbml'])
    producibility_targets_file = os.path.join(*[respath, 'producibility_targets.json'])

    # open and load all result files
    with open(iscope_file, 'r') as json_idata:
        iscope = json.load(json_idata)
    with open(rev_iscope_file, 'r') as json_idata:
        rev_iscope = json.load(json_idata)
    with open(seeds_iscope_file, 'r') as json_idata:
        seeds_iscope = json.load(json_idata)
    with open(cscope_file, 'r') as json_cdata:
        cscope = json.load(json_cdata)
    with open(rev_cscope_file, 'r') as json_cdata:
        rev_cscope = json.load(json_cdata)
    with open(contrib_microbes_file, 'r') as json_cdata:
        contrib_microbes = json.load(json_cdata)
    with open(addedvalue_file, 'r') as json_cdata:
        addedvalue = json.load(json_cdata)
    with open(mincom_file, 'r') as json_cdata:
        mincom = json.load(json_cdata)
    with open(producibility_targets_file, 'r') as json_data:
        producibility_targets = json.load(json_data)

    # ISCOPE ANALYSIS
    # ensure there is the right number of computed indiv scopes
    assert len(iscope) == NUMBER_BACT
    # ensure the union and intersection are ok
    iscope_set = {}
    for elem in iscope:
        iscope_set[elem] = set(iscope[elem])
    # union of iscopes
    union_iscope = set.union(*list(iscope_set.values()))
    assert len(union_iscope) == SIZE_UNION
    intersection_iscope = set.intersection(*list(iscope_set.values()))
    # intersection of iscopes
    assert len(intersection_iscope) == SIZE_INTERSECTION
    # scope content
    for bact in iscope:
        assert set(iscope[bact]) == set(ISCOPE[bact])
    # reverse iscope
    for compound in rev_iscope:
        assert set(rev_iscope[compound]) == set(REV_ISCOPE[compound])
    # seeds in iscope
    for category in PRODUCED_SEEDS_ISCOPE:
        for bact in PRODUCED_SEEDS_ISCOPE[category]:
            assert set(seeds_iscope[category][bact]) == set(PRODUCED_SEEDS_ISCOPE[category][bact])
    
    # CSCOPE ANALYSIS
    # comscope content
    assert set(cscope['com_scope']) == set(COMSCOPE['com_scope'])
    # reverse cscope
    for compound in rev_cscope:
        assert set(rev_cscope[compound]) == set(REV_CSCOPE[compound])
    # contributions of microbes
    for bact in contrib_microbes:
        for category in contrib_microbes[bact]:
            assert set(contrib_microbes[bact][category]) == set(MICROBE_CONTRIBUTIONS[bact][category])

    # ADDEDVALUE ANALYSIS
    # newly producible compounds
    assert set(addedvalue['addedvalue']) == EXPECTED_TARGETS_ADVAL
    reader = SBMLReader()
    document = reader.readSBML(targets_file)
    new_targets = set([specie.getId() for specie in document.getModel().getListOfSpecies()])
    assert new_targets == PROD_TARGETS.union(UNPROD_TARGETS)

    # MINCOM ANALYSIS
    # ensure the minimal number of bacteria in a minimal community is ok
    assert len(mincom['bacteria']) == MIN_SIZE_COM
    # ensure the bacteria in union are ok
    assert set(mincom['union_bacteria']) == UNION_MINCOM
    # ensure the bacteria in intersection are ok
    assert set(mincom['inter_bacteria']) == INTERSECTION_MINCOM
    # ensure the newly producible targets are ok
    assert set(mincom['producible']) == PROD_TARGETS

    # PRODUCIBILITY ANALYSIS
    for key in PRODUCIBILITY_TARGETS:
        if key in ["unproducible", "producible", "indiv_producible", "mincom_producible", "key_species"]:
            assert set(producibility_targets[key]) == set(PRODUCIBILITY_TARGETS[key])
        else:
            for compound in producibility_targets[key]:
                assert set(producibility_targets[key][compound]) == set(PRODUCIBILITY_TARGETS[key][compound])

    # clean
    # Due to unstable behaviour of os.unlink on Windows, do not delete the file.
    # Refer to: https://github.com/python/cpython/issues/109608
    if sys.platform != 'win32':
        shutil.rmtree(respath)


def test_m2m_metacom_tiny_toy_host():
    """
    Test m2m metacom when called from the API using the tiny toy dataset.
    """
    # RUN THE COMMAND
    inppath = 'metabolic_data/tiny_toy'
    respath = 'tiny_metacom_output'
    networks_path = os.path.join(inppath, 'networks')
    seeds_path = os.path.join(inppath, 'seeds_community.sbml')
    targets_path = os.path.join(inppath, 'targets_community.sbml')
    host_path = os.path.join(inppath, 'host.sbml')

    metacom_analysis(sbml_dir=networks_path, out_dir=respath, seeds=seeds_path, host_mn=host_path, targets_file=targets_path, cpu_number=1)

    mincom_file = os.path.join(*[respath, 'community_analysis', 'mincom.json'])
    with open(mincom_file, 'r') as json_cdata:
        mincom = json.load(json_cdata)

    # MINCOM ANALYSIS
    # ensure the minimal number of bacteria in a minimal community is ok
    assert len(mincom['bacteria']) == MIN_SIZE_COM
    # ensure the bacteria in union are ok
    assert set(mincom['union_bacteria']) == UNION_MINCOM
    # ensure the bacteria in intersection are ok
    assert set(mincom['inter_bacteria']) == INTERSECTION_MINCOM
    # ensure the newly producible targets are ok
    assert set(mincom['producible']) == PROD_TARGETS

    comm_scopes_file = os.path.join(*[respath, 'community_analysis', 'comm_scopes.json'])
    with open(comm_scopes_file, 'r') as json_cdata:
        comm_scopes = json.load(json_cdata)

    # Check host keys in comm_scopes json.
    assert len(comm_scopes['com_scope']) == 13
    assert set(comm_scopes['host_unprodtargets']) == set(['M_F_c', 'M_H_c', 'M_foo_c', 'M_C_c'])
    # Host scope = old behavior, seed in scope.
    assert set(comm_scopes['host_scope']) == set(['M_S1_c', 'M_S2_c'])

    # clean
    # Due to unstable behaviour of os.unlink on Windows, do not delete the file.
    # Refer to: https://github.com/python/cpython/issues/109608
    if sys.platform != 'win32':
        shutil.rmtree(respath)


if __name__ == "__main__":
    test_m2m_metacom_tiny_toy()
