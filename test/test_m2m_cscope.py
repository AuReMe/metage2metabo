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
from libsbml import SBMLReader


EXPECTED_TARGETS = {
        'M_3__45__OCTAPRENYL__45__4__45__HYDROXYBENZOATE_c',
        'M_CPD__45__16607_c', 'M_CPD__45__221_c', 'M_CPD__45__16020_c',
        'M_D__45__LACTOYL__45__COA_c', 'M_CPD__45__7695_c',
        'M_CPD0__45__2331_c', 'M_CARBON__45__MONOXIDE_c', 'M_CPD__45__12306_c',
        'M_CPD0__45__2123_c', 'M_HSO3_c', 'M_CPD__45__12179_c',
        'M_CPD0__45__2338_c', 'M_CPD__45__641_c', 'M_CPD__45__421_c',
        'M_ACETOACETYL__45__COA_c', 'M_FADH2_c', 'M_TREHALOSE__45__6P_c',
        'M_ACETONE_c', 'M_CPD__45__592_c', 'M_DEOXYCYTIDINE_c',
        'M_C55__45__PP__45__GLCNAC__45__MANNACA_c',
        'M_ACETYL__45__ETCETERA__45__GLUCOSAMINYLDIPHOSPHOUND_c',
        'M_HOMO__45__SER_c', 'M_RETINAL_c',
        'M_S__45__3__45__HYDROXYBUTANOYL__45__COA_c', 'M_CPD__45__672_c',
        'M_CPD__45__18529_c',
        'M_5__45__PHOSPHORIBOSYL__45__N__45__FORMYLGLYCINEAMIDINE_c',
        'M_DELTA3__45__ISOPENTENYL__45__PP_c', 'M_CPD__45__507_c',
        'M_K__45__HEXANOYL__45__COA_c', 'M_OCTAPRENYL__45__DIPHOSPHATE_c',
        'M_CMP__45__KDO_c', 'M_DEOXYADENOSINE_c', 'M_OXALYL__45__COA_c',
        'M_C4_c', 'M_METHACRYLYL__45__COA_c', 'M_FARNESYL__45__PP_c',
        'M_PHOSPHORIBOSYL__45__CARBOXY__45__AMINOIMIDAZOLE_c',
        'M_CROTONYL__45__COA_c', 'M_MEVALONATE_c',
        'M_4__45__GUANIDO__45__BUTYRAMIDE_c', 'M_CPD__45__12305_c',
        'M_DEOXYURIDINE_c', 'M_2__45__OCTAPRENYL__45__6__45__METHOXYPHENOL_c',
        'M_CPD__45__650_c', 'M_CPD__45__14017_c', 'M_CPD__45__14021_c',
        'M_5__45__PHOSPHORIBOSYL__45__5__45__AMINOIMIDAZOLE_c',
        'M_INDOLEYL__45__CPD_c', 'M_CPD__45__569_c',
        'M_O__45__PHOSPHO__45__L__45__HOMOSERINE_c', 'M_CPD__45__7000_c',
        'M_CPD__45__568_c', 'M_LACTOYL__45__COA_c', 'M_CPD__45__12303_c',
        'M_CPD__45__12173_c', 'M_FORMYL__45__COA_c', 'M_CPD__45__499_c',
        'M_CPD0__45__2340_c',
        'M_OCTAPRENYL__45__METHYL__45__METHOXY__45__BENZQ_c',
        'M_2__45__OCTAPRENYL__45__6__45__HYDROXYPHENOL_c',
        'M_GERANYL__45__PP_c', 'M_Teichoic__45__P__45__Gro__45__Glc_c',
        'M_CPD__45__4211_c', 'M_CPD__45__9646_c', 'M_CPD__45__12177_c',
        'M_OCTAPRENYL__45__METHOXY__45__BENZOQUINONE_c', 'M_CPD__45__466_c',
        'M_DADP_c', 'M_CPD__45__5802_c', 'M_CPD__45__804_c',
        'M_ACRYLYL__45__COA_c', 'M_UNDECAPRENYL__45__DIPHOSPHATE_c',
        'M_Teichoic__45__P__45__Gro_c',
        'M_2__45__KETO__45__3__45__DEOXY__45__D__45__GLUCARATE_c',
        'M_CPD__45__12307_c', 'M_OH__45__HEXANOYL__45__COA_c',
        'M_CPD0__45__181_c', 'M_CPD__45__822_c',
        'M_Teichoic__45__P__45__Gro__45__Glc_e', 'M_CPD0__45__2339_c',
        'M_CPD__45__12304_c', 'M_CPD__45__448_c', 'M_CPD__45__3462_c',
        'M_ACETYL__45__D__45__GLUCOSAMINYLDIPHOSPHO__45__UNDECAPRE_c',
        'M_2__45__OCTAPRENYLPHENOL_c', 'M_UDP__45__GLUCURONATE_c',
        'M_CPD__45__9852_c', 'M_DAMP_c', 'M_CPD__45__15199_c',
        'M_CPD0__45__2244_c',
        'M_UDP__45__4__45__AMINO__45__4__45__DEOXY__45__L__45__ARABINOSE_c',
        'M_CPD__45__12117_c', 'M_S2O3_c', 'M_CPD__45__12175_c',
        'M_CPD__45__12125_c', 'M_CPD__45__19958_c',
        'M_5__45__P__45__RIBOSYL__45__N__45__FORMYLGLYCINEAMIDE_c',
        'M_CPD__45__12115_c', 'M_N2__45__SUCCINYLORNITHINE_c',
        'M_CPD__45__597_c',
        'M_ALL__45__TRANS__45__HEPTAPRENYL__45__DIPHOSPHATE_c',
        'M_3__45__HYDROXY__45__3__45__METHYL__45__GLUTARYL__45__COA_c',
        'M_TREHALOSE_c', 'M_CPD__45__12309_c',
        'M_CH3__45__MALONATE__45__S__45__ALD_c', 'M_GLC__45__6__45__P_c',
        'M_CPD__45__12308_c', 'M_CPD__45__9853_c',
        'M_5__45__BETA__45__L__45__THREO__45__PENTAPYRANOSYL__45__4__45__ULOSE__45___c',
        'M_REDUCED__45__MENAQUINONE_c', 'M_CPD__45__14378_c',
        'M_2__45__METHYL__45__ACETO__45__ACETYL__45__COA_c',
        'M_CPD0__45__2121_c', 'M_CPD__45__12773_c', 'M_CPD__45__13644_c',
        'M_ALPHA__45__GLC__45__6__45__P_c'
    }

SIZE_CSCOPE = 651


def test_m2m_cscope_call():
    """
    Test m2m addedvalue when called in terminal.
    """
    # RUN THE COMMAND
    inppath = 'metabolic_data'
    respath = 'addedvalue_output'
    toy_bact_tgz = os.path.join(inppath, 'toy_bact.tar.gz')
    toy_bact_path = os.path.join(respath, 'toy_bact')
    seeds_path = os.path.join(inppath, 'seeds_toy.sbml')
    targets_path = os.path.join(inppath, 'targets_toy.sbml')

    if not os.path.exists(respath):
        os.makedirs(respath)
    with tarfile.open(toy_bact_tgz) as tar:
        tar.extractall(path=respath)
    subprocess.call([
        'm2m', 'cscope', '-n', toy_bact_path, '-o',
        respath, '-s', seeds_path,
        '-t', targets_path, '-q'
    ])
    cscope_file = os.path.join(*[respath, 'community_analysis', 'comm_scopes.json'])
    # CSCOPE ANALYSIS
    with open(cscope_file, 'r') as json_cdata:
        d_cscope = json.load(json_cdata)
    assert len(d_cscope['com_scope']) == SIZE_CSCOPE
    assert sorted(d_cscope['com_prodtargets']) == sorted(EXPECTED_TARGETS)
    assert d_cscope['com_unprodtargets'] == []

    # clean
    shutil.rmtree(respath)

if __name__ == "__main__":
    test_m2m_cscope_call()