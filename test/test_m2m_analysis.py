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

KEY_SPECIES = ['GCA_003437665', 'GCA_003437785', 'GCA_003437345',
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
    inppath = 'metabolic_data'
    respath = 'm2m_analysis_output'
    draft_path = os.path.join(respath, 'toy_bact')
    draft_tgz_path = os.path.join(inppath, 'toy_bact.tar.gz')
    seeds_path = os.path.join(inppath, 'seeds_toy.sbml')
    targets_path = os.path.join(inppath, 'targets_toy.sbml')
    json_folder = os.path.join(respath, 'json')
    json_file = os.path.join(json_folder, 'targets_toy.json')
    gml_file_path = os.path.join(*[respath, 'gml', 'targets_toy.gml'])

    if not os.path.exists(respath):
        os.makedirs(respath)
    with tarfile.open(draft_tgz_path) as tar:
        tar.extractall(path=respath)
    subprocess.call([
        'm2m_analysis', 'enum', '-n', draft_path, '-o',
        respath, '-t', targets_path, '-s', seeds_path,
        '-q'])


    # KEY SPECIES ANALYSIS
    with open(json_file, 'r') as json_data:
        d_m2m_analysis = json.load(json_data)
    assert sorted(d_m2m_analysis['union_bacteria']) == sorted(KEY_SPECIES)
    assert sorted(d_m2m_analysis['inter_bacteria']) == sorted(ESSENTIAL_SYMBIONTS)
    assert sorted(set(d_m2m_analysis['union_bacteria'])-set(d_m2m_analysis['inter_bacteria'])) == sorted(ALTERNATIVE_SYMBIONTS)
    expected_solutions = sorted([sorted(sol) for sol in list(ENUM_BACTERIA.values())])
    found_solutions = sorted([sorted(sol) for sol in list(d_m2m_analysis['enum_bacteria'].values())])
    assert found_solutions == expected_solutions

    subprocess.call([
        'm2m_analysis', 'graph', '-j', json_folder,
        '-o', respath, '-t', targets_path, '-q'])

    # Graph analysis
    G = nx.read_gml(gml_file_path)
    assert sorted(G.nodes) == sorted(GML_NODES)
    assert len(G.edges()) == 126
    # clean
    shutil.rmtree(respath)

if __name__ == "__main__":
    test_m2m_analysis_call()