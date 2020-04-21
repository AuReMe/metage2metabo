#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test m2m addedvalue on 17 metabolic networks and a file representing growth medium (seeds).
"""

import os
import shutil
import subprocess
import tarfile
import json
import networkx as nx

KEYSTONE_SPECIES = ['GCA_003437665', 'GCA_003437785', 'GCA_003437345',
                     'GCA_003437325', 'GCA_003437055', 'GCA_003437195',
                     'GCA_003437375', 'GCA_003437905', 'GCA_003437715',
                     'GCA_003437885', 'GCA_003437945', 'GCA_003437815',
                     'GCA_003437295', 'GCA_003437595', 'GCA_003437175',
                     'GCA_003437255', 'GCA_003438055']

ESSENTIAL_SYMBIONTS = ['GCA_003437665', 'GCA_003437905', 'GCA_003437715',
                        'GCA_003437885', 'GCA_003437055', 'GCA_003437815',
                        'GCA_003437295', 'GCA_003437595', 'GCA_003437375',
                        'GCA_003437255', 'GCA_003438055', 'GCA_003437195']

ALTERNATIVE_SYMBIONTS = ['GCA_003437345', 'GCA_003437175', 'GCA_003437325',
                            'GCA_003437945', 'GCA_003437785']

ENUM_BACTERIA = {'1': ['GCA_003437905', 'GCA_003437375', 'GCA_003437255',
                        'GCA_003437665', 'GCA_003437195', 'GCA_003437175',
                        'GCA_003437885', 'GCA_003437715', 'GCA_003437815',
                        'GCA_003437295', 'GCA_003438055', 'GCA_003437595',
                        'GCA_003437055'],
                 '2': ['GCA_003437905', 'GCA_003437375', 'GCA_003437255',
                        'GCA_003437345', 'GCA_003437665', 'GCA_003437195',
                        'GCA_003437885', 'GCA_003437715', 'GCA_003437815',
                        'GCA_003437295', 'GCA_003438055', 'GCA_003437595',
                        'GCA_003437055'],
                 '3': ['GCA_003437905', 'GCA_003437375', 'GCA_003437255',
                        'GCA_003437665', 'GCA_003437195', 'GCA_003437885',
                        'GCA_003437715', 'GCA_003437815', 'GCA_003437785',
                        'GCA_003437295', 'GCA_003438055', 'GCA_003437595',
                        'GCA_003437055'],
                 '4': ['GCA_003437905', 'GCA_003437055', 'GCA_003437375',
                        'GCA_003437255', 'GCA_003437665', 'GCA_003437195',
                        'GCA_003437885', 'GCA_003437715', 'GCA_003437815',
                        'GCA_003437295', 'GCA_003438055', 'GCA_003437595',
                        'GCA_003437325'],
                 '5': ['GCA_003437905', 'GCA_003437945', 'GCA_003437375',
                        'GCA_003437255', 'GCA_003437665', 'GCA_003437195',
                        'GCA_003437885', 'GCA_003437715', 'GCA_003437815',
                        'GCA_003437295', 'GCA_003438055', 'GCA_003437595',
                        'GCA_003437055']
                }

GML_NODES = ['GCA_003437595', 'GCA_003437055', 'GCA_003437665', 
             'GCA_003437715', 'GCA_003437815', 'GCA_003437885',
             'GCA_003437905', 'GCA_003437255', 'GCA_003437375',
             'GCA_003437295', 'GCA_003438055', 'GCA_003437175',
             'GCA_003437195', 'GCA_003437345', 'GCA_003437785',
             'GCA_003437325', 'GCA_003437945']

def test_m2m_analysis_call():
    """
    Test m2m analysis when called in terminal.
    """
    # RUN THE COMMAND
    inppath = 'metabolic_data/'
    respath = 'm2m_analysis_output/'
    if not os.path.exists(respath):
        os.makedirs(respath)
    with tarfile.open(inppath + 'toy_bact.tar.gz') as tar:
        tar.extractall(path=respath)
    subprocess.call([
        'm2m_analysis', 'enum', '-n', respath + '/toy_bact', '-o',
        respath, '-t', inppath + '/targets_toy.sbml', '-s', inppath + '/seeds_toy.sbml',
        '-q'])

    json_file = respath +'json' + '/' + 'targets_toy.json'

    # KEYSTONE SPECIES ANALYSIS
    with open(json_file, 'r') as json_data:
        d_m2m_analysis = json.load(json_data)
    assert sorted(d_m2m_analysis['union_bacteria']) == sorted(KEYSTONE_SPECIES)
    assert sorted(d_m2m_analysis['inter_bacteria']) == sorted(ESSENTIAL_SYMBIONTS)
    assert sorted(set(d_m2m_analysis['union_bacteria'])-set(d_m2m_analysis['inter_bacteria'])) == sorted(ALTERNATIVE_SYMBIONTS)
    expected_solutions = sorted([sorted(sol) for sol in list(ENUM_BACTERIA.values())])
    found_solutions = sorted([sorted(sol) for sol in list(d_m2m_analysis['enum_bacteria'].values())])
    assert found_solutions == expected_solutions

    subprocess.call([
        'm2m_analysis', 'graph', '-j', respath + '/json',
        '-o', respath, '-t', inppath + '/targets_toy.sbml', '-q'])
    gml_file = respath +'gml' + '/' + 'targets_toy.gml'

    # Graph analysis
    G = nx.read_gml(gml_file)
    assert sorted(G.nodes) == sorted(GML_NODES)
    assert len(G.edges()) == 126
    # clean
    shutil.rmtree(respath)

if __name__ == "__main__":
    test_m2m_analysis_call()