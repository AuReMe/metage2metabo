#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Test m2m iscope on 17 metabolic networks and a file representing growth medium (seeds).
"""

import csv
import os
import shutil
import subprocess
import tarfile
import json

EXPECTED_PRODUCED_COMPOUNDS = {
    'M_10__45__FORMYL__45__THF_c': 17,
    'M_1__45__AMINO__45__PROPAN__45__2__45__ONE__45__3__45__PHOSPHATE_c': 13,
    'M_1__45__KESTOTRIOSE_c': 1,
    'M_1__45__KETO__45__2__45__METHYLVALERATE_c': 13,
    'M_1__45__L__45__MYO__45__INOSITOL__45__1__45__P_c': 7,
    'M_2C__45__METH__45__D__45__ERYTHRITOL__45__CYCLODIPHOSPHATE_c': 9,
    'M_2K__45__4CH3__45__PENTANOATE_c': 3,
    'M_2__45__3__45__DIHYDROXYBENZOATE_c': 2,
    'M_2__45__5__45__TRIPHOSPHORIBOSYL__45__3__45__DEPHOSPHO__45___c': 1,
    'M_2__45__ACETO__45__2__45__HYDROXY__45__BUTYRATE_c': 13,
    'M_2__45__ACETO__45__LACTATE_c': 15,
    'M_2__45__ALPHA__45__HYDROXYETHYL__45__THPP_c': 4,
    'M_2__45__AMINOACRYLATE_c': 16,
    'M_2__45__C__45__METHYL__45__D__45__ERYTHRITOL__45__4__45__PHOSPHATE_c': 15,
    'M_2__45__D__45__THREO__45__HYDROXY__45__3__45__CARBOXY__45__ISOCAPROATE_c': 2,
    'M_2__45__HYDROXY__45__3__45__KETO__45__5__45__METHYLTHIO__45__1__45__PHOSPHOP_c': 2,
    'M_2__45__KETOGLUTARATE_c': 16,
    'M_2__45__KETO__45__3__45__DEOXY__45__6__45__P__45__GLUCONATE_c': 3,
    'M_2__45__KETO__45__3__45__METHYL__45__VALERATE_c': 15,
    'M_2__45__KETO__45__GLUTARAMATE_c': 2,
    'M_2__45__KETO__45__ISOVALERATE_c': 15,
    'M_2__45__METHYL__45__3__45__PHYTYL__45__14__45__NAPHTHOQUINONE_c': 17,
    'M_2__45__METHYL__45__BUTYRYL__45__COA_c': 1,
    'M_2__45__OXOBUTANOATE_c': 13,
    'M_2__45__PG_c': 17,
    'M_2__45__PHOSPHO__45__4__45__CYTIDINE__45__5__45__DIPHOSPHO__45__2__45__C__45__MET_c': 9,
    'M_3OH__45__4P__45__OH__45__ALPHA__45__KETOBUTYRATE_c': 13,
    'M_3S__45__CITRYL__45__COA_c': 1,
    'M_3__45__4__45__DIHYDROXYBENZOATE_c': 1,
    'M_3__45__CARBOXY__45__3__45__HYDROXY__45__ISOCAPROATE_c': 2,
    'M_3__45__DEHYDRO__45__SHIKIMATE_c': 15,
    'M_3__45__DEOXY__45__D__45__ARABINO__45__HEPTULOSONATE__45__7__45__P_c': 15,
    'M_3__45__ENOLPYRUVYL__45__SHIKIMATE__45__5P_c': 15,
    'M_3__45__HYDROXY__45__PROPIONATE_c': 2,
    'M_3__45__KETOBUTYRATE_c': 1,
    'M_3__45__KETO__45__ADIPATE_c': 1,
    'M_3__45__KETO__45__L__45__GULONATE_c': 3,
    'M_3__45__MERCAPTO__45__PYRUVATE_c': 2,
    'M_3__45__OXOADIPATE__45__ENOL__45__LACTONE_c': 1,
    'M_3__45__P__45__HYDROXYPYRUVATE_c': 14,
    'M_3__45__P__45__SERINE_c': 14,
    'M_3__45__SULFINOALANINE_c': 3,
    'M_3__45__UREIDO__45__PROPIONATE_c': 2,
    'M_4__45__AMINO__45__4__45__DEOXYCHORISMATE_c': 13,
    'M_4__45__AMINO__45__BUTYRALDEHYDE_c': 2,
    'M_4__45__AMINO__45__BUTYRATE_c': 5,
    'M_4__45__CYTIDINE__45__5__45__DIPHOSPHO__45__2__45__C_c': 9,
    'M_4__45__FUMARYL__45__ACETOACETATE_c': 1,
    'M_4__45__IMIDAZOLONE__45__5__45__PROPIONATE_c': 4,
    'M_4__45__MALEYL__45__ACETOACETATE_c': 1,
    'M_4__45__PHOSPHONOOXY__45__THREONINE_c': 13,
    'M_4__45__P__45__PANTOTHENATE_c': 8,
    'M_4__45__hydroxybenzoate_c': 1,
    'M_5__45__HYDROXYISOURATE_c': 2,
    'M_5__45__HYDROXY__45__CTP_c': 4,
    'M_5__45__METHYLTHIOADENOSINE_c': 5,
    'M_5__45__METHYL__45__THF_c': 17,
    'M_5__45__OXOPROLINE_c': 15,
    'M_5__45__PHOSPHO__45__RIBOSYL__45__GLYCINEAMIDE_c': 15,
    'M_5__45__P__45__BETA__45__D__45__RIBOSYL__45__AMINE_c': 15,
    'M_6__45__KESTOSE_c': 1,
    'M_7__45__8__45__DIHYDROPTEROATE_c': 5,
    'M_7__45__AMINOMETHYL__45__7__45__DEAZAGUANINE_c': 5,
    'M_7__45__CYANO__45__7__45__DEAZAGUANINE_c': 5,
    'M_ACETALD_c': 10,
    'M_ACETYLSERINE_c': 3,
    'M_ACETYL__45__COA_c': 3,
    'M_ACETYL__45__GLU_c': 1,
    'M_ACETYL__45__P_c': 3,
    'M_ACET_c': 6,
    'M_ADENINE_c': 9,
    'M_ADENOSINE_c': 16,
    'M_ADENOSYLCOBALAMIN_c': 17,
    'M_ADENOSYL__45__HOMO__45__CYS_c': 3,
    'M_ADENOSYL__45__P4_c': 2,
    'M_ADENYLOSUCC_c': 10,
    'M_ADP__45__D__45__GLUCOSE_c': 13,
    'M_ADP__45__D__45__GLYCERO__45__D__45__MANNO__45__HEPTOSE_c': 2,
    'M_ADP__45__L__45__GLYCERO__45__D__45__MANNO__45__HEPTOSE_c': 2,
    'M_ADP_c': 17,
    'M_AGMATHINE_c': 6,
    'M_AICAR_c': 14,
    'M_ALLANTOATE_c': 2,
    'M_ALPHA__45__GLUCOSE_c': 17,
    'M_ALPHA__45__TOCOPHEROL_c': 17,
    'M_AMINO__45__ACETONE_c': 6,
    'M_AMINO__45__HYDROXYMETHYL__45__METHYLPYRIMIDINE__45__PP_c': 4,
    'M_AMINO__45__HYDROXYMETHYL__45__METHYL__45__PYR__45__P_c': 4,
    'M_AMINO__45__OH__45__HYDROXYMETHYL__45__DIHYDROPTERIDINE_c': 6,
    'M_AMINO__45__OXOBUT_c': 6,
    'M_AMINO__45__RIBOSYLAMINO__45__1H__45__3H__45__PYR__45__DIONE_c': 6,
    'M_AMMONIA_c': 17,
    'M_AMMONIUM_c': 17,
    'M_AMP_c': 17,
    'M_ANTHRANILATE_c': 13,
    'M_ARABINOSE__45__5P_c': 2,
    'M_ARACHIDIC_ACID_c': 17,
    'M_ARACHIDONIC_ACID_c': 17,
    'M_ARG_c': 17,
    'M_ASCORBATE_c': 17,
    'M_ASN_c': 9,
    'M_ATP_c': 17,
    'M_BETA__45__D__45__FRUCTOSE_c': 17,
    'M_BIFURCOSE_c': 1,
    'M_BIOTIN_c': 17,
    'M_BIO__45__5__45__AMP_c': 14,
    'M_BUTYRIC_ACID_c': 17,
    'M_BUTYRYL__45__COA_c': 2,
    'M_BUTYRYL__45__P_c': 6,
    'M_B__45__ALANINE_c': 7,
    'M_C3_c': 1,
    'M_CADAVERINE_c': 2,
    'M_CAMP_c': 4,
    'M_CARBAMATE_c': 13,
    'M_CARBAMOYL__45__P_c': 13,
    'M_CARBAMYUL__45__L__45__ASPARTATE_c': 13,
    'M_CARBON__45__DIOXIDE_c': 17,
    'M_CARBOXYPHENYLAMINO__45__DEOXYRIBULOSE__45__P_c': 13,
    'M_CA__43__2_c': 17,
    'M_CDP__45__D__45__GLUCOSE_c': 3,
    'M_CDP_c': 10,
    'M_CELLOBIOSE_c': 1,
    'M_CELLULOSE_c': 17,
    'M_CGMP_c': 1,
    'M_CH33ADO_c': 12,
    'M_CHOLESTEROL_c': 17,
    'M_CHORISMATE_c': 15,
    'M_CIS__45__ACONITATE_c': 10,
    'M_CIT_c': 12,
    'M_CL__45___c': 17,
    'M_CMP_c': 10,
    'M_CO3_c': 13,
    'M_CO__45__A_c': 3,
    'M_CPD0__45__1108_c': 10,
    'M_CPD0__45__1456_c': 3,
    'M_CPD0__45__1699_c': 5,
    'M_CPD0__45__1905_c': 9,
    'M_CPD0__45__2015_c': 1,
    'M_CPD0__45__2101_c': 2,
    'M_CPD0__45__2208_c': 17,
    'M_CPD0__45__2298_c': 2,
    'M_CPD0__45__2461_c': 9,
    'M_CPD0__45__2467_c': 2,
    'M_CPD0__45__2468_c': 2,
    'M_CPD0__45__2472_c': 17,
    'M_CPD0__45__2474_c': 10,
    'M_CPD0__45__2483_c': 2,
    'M_CPD1F__45__129_c': 17,
    'M_CPD__45__10244_c': 17,
    'M_CPD__45__10267_c': 2,
    'M_CPD__45__10330_c': 10,
    'M_CPD__45__10353_c': 9,
    'M_CPD__45__10551_c': 2,
    'M_CPD__45__1063_c': 3,
    'M_CPD__45__1067_c': 1,
    'M_CPD__45__10750_c': 1,
    'M_CPD__45__10751_c': 1,
    'M_CPD__45__10752_c': 1,
    'M_CPD__45__10753_c': 1,
    'M_CPD__45__10754_c': 1,
    'M_CPD__45__10755_c': 1,
    'M_CPD__45__10773_c': 3,
    'M_CPD__45__10774_c': 3,
    'M_CPD__45__10775_c': 3,
    'M_CPD__45__10776_c': 3,
    'M_CPD__45__1086_c': 6,
    'M_CPD__45__108_c': 12,
    'M_CPD__45__1091_c': 2,
    'M_CPD__45__11281_c': 2,
    'M_CPD__45__1133_c': 2,
    'M_CPD__45__11770_c': 2,
    'M_CPD__45__11855_c': 2,
    'M_CPD__45__12189_c': 17,
    'M_CPD__45__12258_c': 1,
    'M_CPD__45__12279_c': 12,
    'M_CPD__45__12365_c': 3,
    'M_CPD__45__12366_c': 9,
    'M_CPD__45__12367_c': 3,
    'M_CPD__45__12377_c': 9,
    'M_CPD__45__12427_c': 3,
    'M_CPD__45__12575_c': 8,
    'M_CPD__45__12601_c': 4,
    'M_CPD__45__12653_c': 17,
    'M_CPD__45__12826_c': 17,
    'M_CPD__45__13043_c': 5,
    'M_CPD__45__13118_c': 1,
    'M_CPD__45__13357_c': 14,
    'M_CPD__45__13469_c': 17,
    'M_CPD__45__13489_c': 5,
    'M_CPD__45__13792_c': 17,
    'M_CPD__45__13847_c': 1,
    'M_CPD__45__13851_c': 9,
    'M_CPD__45__13907_c': 3,
    'M_CPD__45__13910_c': 3,
    'M_CPD__45__13912_c': 3,
    'M_CPD__45__13913_c': 3,
    'M_CPD__45__13914_c': 3,
    'M_CPD__45__14292_c': 17,
    'M_CPD__45__14443_c': 16,
    'M_CPD__45__14553_c': 8,
    'M_CPD__45__14808_c': 2,
    'M_CPD__45__15015_c': 3,
    'M_CPD__45__15016_c': 3,
    'M_CPD__45__15056_c': 13,
    'M_CPD__45__15127_c': 2,
    'M_CPD__45__15158_c': 2,
    'M_CPD__45__15317_c': 10,
    'M_CPD__45__15358_c': 2,
    'M_CPD__45__15373_c': 4,
    'M_CPD__45__15377_c': 4,
    'M_CPD__45__15382_c': 7,
    'M_CPD__45__15403_c': 2,
    'M_CPD__45__15435_c': 13,
    'M_CPD__45__15590_c': 15,
    'M_CPD__45__155_c': 1,
    'M_CPD__45__15700_c': 3,
    'M_CPD__45__15709_c': 17,
    'M_CPD__45__15818_c': 10,
    'M_CPD__45__15972_c': 17,
    'M_CPD__45__15975_c': 3,
    'M_CPD__45__15979_c': 2,
    'M_CPD__45__159_c': 1,
    'M_CPD__45__16013_c': 13,
    'M_CPD__45__16015_c': 16,
    'M_CPD__45__16606_c': 1,
    'M_CPD__45__16876_c': 2,
    'M_CPD__45__17188_c': 17,
    'M_CPD__45__17322_c': 17,
    'M_CPD__45__18085_c': 17,
    'M_CPD__45__18118_c': 2,
    'M_CPD__45__18238_c': 13,
    'M_CPD__45__187_c': 2,
    'M_CPD__45__18_c': 2,
    'M_CPD__45__19306_c': 2,
    'M_CPD__45__195_c': 17,
    'M_CPD__45__196_c': 2,
    'M_CPD__45__19753_c': 2,
    'M_CPD__45__20826_c': 8,
    'M_CPD__45__209_c': 1,
    'M_CPD__45__227_c': 1,
    'M_CPD__45__233_c': 1,
    'M_CPD__45__2343_c': 3,
    'M_CPD__45__235_c': 1,
    'M_CPD__45__245_c': 1,
    'M_CPD__45__2961_c': 15,
    'M_CPD__45__302_c': 2,
    'M_CPD__45__316_c': 9,
    'M_CPD__45__317_c': 2,
    'M_CPD__45__318_c': 3,
    'M_CPD__45__334_c': 3,
    'M_CPD__45__335_c': 1,
    'M_CPD__45__3617_c': 17,
    'M_CPD__45__365_c': 2,
    'M_CPD__45__375_c': 7,
    'M_CPD__45__380_c': 2,
    'M_CPD__45__387_c': 17,
    'M_CPD__45__389_c': 2,
    'M_CPD__45__444_c': 3,
    'M_CPD__45__469_c': 1,
    'M_CPD__45__470_c': 2,
    'M_CPD__45__479_c': 2,
    'M_CPD__45__551_c': 2,
    'M_CPD__45__560_c': 5,
    'M_CPD__45__564_c': 3,
    'M_CPD__45__578_c': 1,
    'M_CPD__45__5821_c': 2,
    'M_CPD__45__602_c': 6,
    'M_CPD__45__606_c': 3,
    'M_CPD__45__6124_c': 2,
    'M_CPD__45__622_c': 1,
    'M_CPD__45__62_c': 2,
    'M_CPD__45__645_c': 2,
    'M_CPD__45__653_c': 17,
    'M_CPD__45__658_c': 3,
    'M_CPD__45__667_c': 1,
    'M_CPD__45__6972_c': 2,
    'M_CPD__45__69_c': 13,
    'M_CPD__45__7100_c': 2,
    'M_CPD__45__7224_c': 2,
    'M_CPD__45__7671_c': 10,
    'M_CPD__45__7737_c': 1,
    'M_CPD__45__7830_c': 17,
    'M_CPD__45__7836_c': 17,
    'M_CPD__45__8050_c': 2,
    'M_CPD__45__8052_c': 2,
    'M_CPD__45__8259_c': 16,
    'M_CPD__45__8268_c': 1,
    'M_CPD__45__827_c': 2,
    'M_CPD__45__8462_c': 17,
    'M_CPD__45__85_c': 2,
    'M_CPD__45__8999_c': 2,
    'M_CPD__45__9000_c': 2,
    'M_CPD__45__9245_c': 17,
    'M_CPD__45__9451_c': 2,
    'M_CPD__45__9550_c': 1,
    'M_CPD__45__9923_c': 4,
    'M_CPD__45__9924_c': 4,
    'M_CPD__45__9925_c': 2,
    'M_CTP_c': 10,
    'M_CU__43__2_c': 17,
    'M_CYS__45__GLY_c': 9,
    'M_CYS_c': 17,
    'M_CYTIDINE_c': 10,
    'M_C__45__DI__45__GMP_c': 2,
    'M_Cellodextrins_c': 4,
    'M_DATP_c': 9,
    'M_DCDP_c': 2,
    'M_DCMP_c': 2,
    'M_DCTP_c': 4,
    'M_DEAMIDO__45__NAD_c': 16,
    'M_DEHYDROQUINATE_c': 15,
    'M_DEOXYGUANOSINE_c': 8,
    'M_DEOXYINOSINE_c': 8,
    'M_DEOXYXYLULOSE__45__5P_c': 16,
    'M_DEOXY__45__D__45__RIBOSE__45__1__45__PHOSPHATE_c': 8,
    'M_DEOXY__45__RIBOSE__45__5P_c': 8,
    'M_DEPHOSPHO__45__COA_c': 3,
    'M_DGDP_c': 1,
    'M_DGMP_c': 1,
    'M_DGTP_c': 9,
    'M_DIAMINO__45__OH__45__PHOSPHORIBOSYLAMINO__45__PYR_c': 7,
    'M_DIHYDROFOLATE_c': 5,
    'M_DIHYDROMONAPTERIN__45__TRIPHOSPHATE_c': 2,
    'M_DIHYDRONEOPTERIN__45__P3_c': 8,
    'M_DIHYDRONEOPTERIN__45__P_c': 6,
    'M_DIHYDROPTERIN__45__CH2OH__45__PP_c': 6,
    'M_DIHYDROXYACETONE_c': 9,
    'M_DIHYDROXYNAPHTHOATE_c': 2,
    'M_DIHYDROXYPENTANEDIONE_c': 3,
    'M_DIHYDROXY__45__ACETONE__45__PHOSPHATE_c': 17,
    'M_DIHYDROXY__45__BUTANONE__45__P_c': 7,
    'M_DIHYDRO__45__DIOH__45__BENZOATE_c': 2,
    'M_DIHYDRO__45__NEO__45__PTERIN_c': 6,
    'M_DIMETHYL__45__D__45__RIBITYL__45__LUMAZINE_c': 6,
    'M_DI__45__H__45__OROTATE_c': 13,
    'M_DI__45__H__45__URACIL_c': 2,
    'M_DOCOSANOATE_c': 17,
    'M_DODECANOATE_c': 17,
    'M_DPG_c': 17,
    'M_DUMP_c': 4,
    'M_DUTP_c': 4,
    'M_D__45__6__45__P__45__GLUCONO__45__DELTA__45__LACTONE_c': 13,
    'M_D__45__ALANINE_c': 17,
    'M_D__45__ALA__45__D__45__ALA_c': 17,
    'M_D__45__ALPHABETA__45__D__45__HEPTOSE__45__7__45__PHOSPHATE_c': 3,
    'M_D__45__BETA__45__D__45__HEPTOSE__45__17__45__DIPHOSPHATE_c': 3,
    'M_D__45__BETA__45__D__45__HEPTOSE__45__1__45__P_c': 2,
    'M_D__45__ERYTHRO__45__IMIDAZOLE__45__GLYCEROL__45__P_c': 14,
    'M_D__45__GLT_c': 17,
    'M_D__45__GLUCOSAMINE__45__6__45__P_c': 17,
    'M_D__45__LACTATE_c': 10,
    'M_D__45__METHYL__45__MALONYL__45__COA_c': 2,
    'M_D__45__PROLINE_c': 1,
    'M_D__45__RIBULOSE__45__15__45__P2_c': 2,
    'M_D__45__RIBULOSE__45__1__45__P_c': 3,
    'M_D__45__RIBULOSE_c': 5,
    'M_D__45__Ribofuranose_c': 10,
    'M_D__45__Ribopyranose_c': 10,
    'M_D__45__SEDOHEPTULOSE__45__7__45__P_c': 15,
    'M_D__45__XYLULOSE_c': 3,
    'M_D__45__Xylopyranose_c': 4,
    'M_D__45__arabinofuranose_c': 3,
    'M_D__45__arabinopyranose_c': 3,
    'M_D__45__galactopyranose_c': 17,
    'M_D__45__glucopyranose__45__6__45__phosphate_c': 17,
    'M_D__45__mannopyranose_c': 4,
    'M_ENOL__45__PHENYLPYRUVATE_c': 4,
    'M_ENTEROBACTIN_c': 2,
    'M_ERYTHRONATE__45__4P_c': 13,
    'M_ERYTHROSE__45__4P_c': 15,
    'M_ETHYLENE__45__CMPD_c': 2,
    'M_ETOH_c': 17,
    'M_FAD_c': 16,
    'M_FE__43__2_c': 17,
    'M_FE__43__3_c': 17,
    'M_FMNH2_c': 6,
    'M_FMN_c': 16,
    'M_FORMAMIDE_c': 2,
    'M_FORMATE_c': 9,
    'M_FORMYL__45__ISOGLUTAMINE_c': 2,
    'M_FRUCTOSE__45__16__45__DIPHOSPHATE_c': 17,
    'M_FRUCTOSE__45__6P_c': 17,
    'M_FUM_c': 17,
    'M_G3P_c': 17,
    'M_GAP_c': 17,
    'M_GDP__45__4__45__DEHYDRO__45__6__45__DEOXY__45__D__45__MANNOSE_c': 1,
    'M_GDP__45__4__45__DEHYDRO__45__6__45__L__45__DEOXYGALACTOSE_c': 1,
    'M_GDP__45__MANNOSE_c': 2,
    'M_GDP__45__TP_c': 10,
    'M_GDP_c': 10,
    'M_GLC__45__1__45__P_c': 17,
    'M_GLC_c': 17,
    'M_GLN_c': 17,
    'M_GLT_c': 17,
    'M_GLUCONATE_c': 2,
    'M_GLUCOSAMINE__45__1P_c': 17,
    'M_GLUTATHIONE_c': 9,
    'M_GLUTATHIONYLSPERMIDINE_c': 2,
    'M_GLYCERATE_c': 1,
    'M_GLYCEROL__45__3P_c': 10,
    'M_GLYCOLALDEHYDE_c': 7,
    'M_GLYCOLLATE_c': 4,
    'M_GLYCOL_c': 2,
    'M_GLYOX_c': 13,
    'M_GLY_c': 17,
    'M_GMP_c': 10,
    'M_GTP_c': 10,
    'M_GUANINE_c': 10,
    'M_GUANOSINE__45__5DP__45__3DP_c': 10,
    'M_GUANOSINE_c': 10,
    'M_Glucopyranose_c': 12,
    'M_Glucose_c': 4,
    'M_H2CO3_c': 13,
    'M_HCO3_c': 13,
    'M_HISTAMINE_c': 1,
    'M_HISTIDINAL_c': 14,
    'M_HISTIDINOL_c': 14,
    'M_HIS_c': 17,
    'M_HMP_c': 4,
    'M_HOMOGENTISATE_c': 1,
    'M_HOMO__45__CYS_c': 9,
    'M_HS_c': 9,
    'M_HYDROGEN__45__MOLECULE_c': 3,
    'M_HYDROGEN__45__PEROXIDE_c': 9,
    'M_HYPOXANTHINE_c': 10,
    'M_IDP_c': 2,
    'M_ILE_c': 17,
    'M_IMIDAZOLE__45__ACETOL__45__P_c': 14,
    'M_IMINOASPARTATE_c': 16,
    'M_IMP_c': 10,
    'M_INDOLE_ACETALDEHYDE_c': 2,
    'M_INDOLE_ACETATE_AUXIN_c': 2,
    'M_INDOLE_PYRUVATE_c': 2,
    'M_INDOLE__45__3__45__GLYCEROL__45__P_c': 13,
    'M_INDOLE_c': 13,
    'M_INOSINE_c': 10,
    'M_ISOBUTYRYL__45__COA_c': 1,
    'M_ISOCHORISMATE_c': 4,
    'M_ISOGLUTAMINE_c': 2,
    'M_ISOVALERYL__45__COA_c': 1,
    'M_ITP_c': 7,
    'M_KDO__45__8P_c': 2,
    'M_KDO_c': 2,
    'M_K__43___c': 17,
    'M_LAUROYLCOA__45__CPD_c': 3,
    'M_LEU_c': 17,
    'M_LINOLEIC_ACID_c': 17,
    'M_LINOLENIC_ACID_c': 17,
    'M_LINOLENOYL__45__COA_c': 2,
    'M_LYS_c': 17,
    'M_L__45__1__45__LYSOPHOSPHATIDATE_c': 1,
    'M_L__45__ALPHA__45__ALANINE_c': 17,
    'M_L__45__ARGININO__45__SUCCINATE_c': 13,
    'M_L__45__ASPARTATE__45__SEMIALDEHYDE_c': 17,
    'M_L__45__ASPARTATE_c': 17,
    'M_L__45__BETA__45__ASPARTYL__45__P_c': 17,
    'M_L__45__CITRULLINE_c': 13,
    'M_L__45__CYSTATHIONINE_c': 9,
    'M_L__45__CYSTEATE_c': 3,
    'M_L__45__DEHYDRO__45__ASCORBATE_c': 3,
    'M_L__45__DELTA1__45__PYRROLINE_5__45__CARBOXYLATE_c': 16,
    'M_L__45__DI__45__GMP_c': 2,
    'M_L__45__GAMMA__45__GLUTAMYLCYSTEINE_c': 9,
    'M_L__45__GLUTAMATE_GAMMA__45__SEMIALDEHYDE_c': 16,
    'M_L__45__GLUTAMATE__45__5__45__P_c': 15,
    'M_L__45__GLYCERALDEHYDE__45__3__45__PHOSPHATE_c': 7,
    'M_L__45__HISTIDINOL__45__P_c': 14,
    'M_L__45__LACTATE_c': 15,
    'M_L__45__ORNITHINE_c': 17,
    'M_L__45__RIBULOSE__45__5__45__P_c': 13,
    'M_L__45__THREO__45__3__45__METHYL__45__ASPARTATE_c': 1,
    'M_L__45__XYLULOSE__45__5__45__P_c': 8,
    'M_Large__45__branched__45__glucans_c': 4,
    'M_Long__45__linear__45__glucans_c': 1,
    'M_MALONATE__45__S__45__ALD_c': 2,
    'M_MALONYL__45__COA_c': 3,
    'M_MALTOSE_c': 17,
    'M_MAL_c': 10,
    'M_MANNITOL_c': 17,
    'M_MANNOSE__45__1P_c': 2,
    'M_MANNOSE_c': 1,
    'M_MESACONATE_c': 1,
    'M_METHYL__45__GLYOXAL_c': 9,
    'M_METHYL__45__MALONYL__45__COA_c': 1,
    'M_MET_c': 17,
    'M_MG__43__2_c': 17,
    'M_MN__43__2_c': 17,
    'M_MYO__45__INOSITOL_c': 7,
    'M_N2__45__SUCCINYLGLUTAMATE_c': 1,
    'M_NADH_c': 17,
    'M_NADPH_c': 17,
    'M_NADP_c': 17,
    'M_NAD_c': 17,
    'M_NA__43___c': 17,
    'M_NA__43___e': 3,
    'M_NEOKESTOSE_c': 1,
    'M_NIACINAMIDE_c': 17,
    'M_NIACINE_c': 17,
    'M_NICOTINAMIDE_NUCLEOTIDE_c': 17,
    'M_NICOTINAMIDE_RIBOSE_c': 17,
    'M_NICOTINATE_NUCLEOTIDE_c': 16,
    'M_NMNH_c': 13,
    'M_NYSTOSE_c': 1,
    'M_N__45__23__45__DIHYDROXYBENZOYL__45__L__45__SERINE_c': 2,
    'M_N__45__5__45__PHOSPHORIBOSYL__45__ANTHRANILATE_c': 13,
    'M_N__45__ACETYL__45__D__45__GLUCOSAMINE__45__1__45__P_c': 3,
    'M_N__45__ACETYL__45__GLUTAMYL__45__P_c': 1,
    'M_N__45__ALPHA__45__ACETYLORNITHINE_c': 1,
    'M_N__45__FORMIMINO__45__L__45__GLUTAMATE_c': 4,
    'M_OH__45__PYR_c': 1,
    'M_OH_c': 9,
    'M_OLEATE__45__CPD_c': 17,
    'M_OLEOYL__45__COA_c': 2,
    'M_OROTATE_c': 11,
    'M_OROTIDINE__45__5__45__PHOSPHATE_c': 10,
    'M_OXALACETIC_ACID_c': 17,
    'M_OXALATE_c': 2,
    'M_OXALO__45__SUCCINATE_c': 13,
    'M_OXAMATE_c': 2,
    'M_OXIDIZED__45__GLUTATHIONE_c': 3,
    'M_OXYGEN__45__MOLECULE_c': 9,
    'M_O__45__SUCCINYLBENZOATE_c': 4,
    'M_O__45__SUCCINYL__45__L__45__HOMOSERINE_c': 9,
    'M_P3I_c': 9,
    'M_PALMITATE_c': 17,
    'M_PALMITYL__45__COA_c': 3,
    'M_PANTETHEINE__45__P_c': 3,
    'M_PANTOTHENATE_c': 17,
    'M_PHENYL__45__PYRUVATE_c': 4,
    'M_PHE_c': 17,
    'M_PHOSPHORIBOSYL__45__AMP_c': 14,
    'M_PHOSPHORIBOSYL__45__ATP_c': 14,
    'M_PHOSPHORIBOSYL__45__FORMAMIDO__45__CARBOXAMIDE_c': 1,
    'M_PHOSPHORIBOSYL__45__FORMIMINO__45__AICAR__45__P_c': 14,
    'M_PHOSPHORIBULOSYL__45__FORMIMINO__45__AICAR__45__P_c': 14,
    'M_PHOSPHO__45__ENOL__45__PYRUVATE_c': 17,
    'M_PPI_c': 17,
    'M_PREPHENATE_c': 10,
    'M_PROPIONATE_c': 1,
    'M_PROPIONYL__45__COA_c': 1,
    'M_PROPIONYL__45__P_c': 1,
    'M_PROTON_c': 17,
    'M_PROTON_e': 17,
    'M_PRO_c': 17,
    'M_PRPP_c': 16,
    'M_PUTRESCINE_c': 6,
    'M_PYRIDOXAL_PHOSPHATE_c': 15,
    'M_PYRIDOXAL_c': 17,
    'M_PYRIDOXAMINE__45__5P_c': 15,
    'M_PYRIDOXAMINE_c': 17,
    'M_PYRIDOXINE__45__5P_c': 15,
    'M_PYRIDOXINE_c': 17,
    'M_PYRUVATE_c': 17,
    'M_P__45__AMINO__45__BENZOATE_c': 13,
    'M_P__45__HYDROXY__45__PHENYLPYRUVATE_c': 3,
    'M_P__45__RIBOSYL__45__4__45__SUCCCARB__45__AMINOIMIDAZOLE_c': 1,
    'M_Pi_c': 17,
    'M_QUINATE_c': 1,
    'M_QUINOLINATE_c': 16,
    'M_RIBOFLAVIN_c': 17,
    'M_RIBOSE__45__1P_c': 10,
    'M_RIBOSE__45__5P_c': 16,
    'M_RIBULOSE__45__5P_c': 16,
    'M_R__45__4__45__PHOSPHOPANTOTHENOYL__45__L__45__CYSTEINE_c': 3,
    'M_R__45____45__ALLANTOIN_c': 2,
    'M_Retinols_c': 17,
    'M_SERYL__45__AMP_c': 2,
    'M_SER_c': 17,
    'M_SHIKIMATE__45__5P_c': 15,
    'M_SHIKIMATE_c': 15,
    'M_SO3_c': 1,
    'M_SORBITOL_c': 17,
    'M_SPERMIDINE_c': 5,
    'M_STEARIC_ACID_c': 17,
    'M_STEAROYL__45__COA_c': 2,
    'M_SUCC__45__S__45__ALD_c': 2,
    'M_SUCROSE_c': 17,
    'M_SUC__45__COA_c': 1,
    'M_SUC_c': 15,
    'M_SUPER__45__OXIDE_c': 8,
    'M_S__45__ADENOSYLMETHIONINAMINE_c': 6,
    'M_S__45__ADENOSYLMETHIONINE_c': 17,
    'M_S__45__ALLANTOIN_c': 2,
    'M_S__45__CITRAMALATE_c': 1,
    'M_S__45__LACTOYL__45__GLUTATHIONE_c': 3,
    'M_Short__45__glucans_c': 1,
    'M_Starch_c': 17,
    'M_TARTRONATE__45__S__45__ALD_c': 1,
    'M_TETRACOSANOATE_c': 17,
    'M_THF_c': 17,
    'M_THIAMINE__45__PYROPHOSPHATE_c': 16,
    'M_THIAMINE__45__P_c': 4,
    'M_THIAMINE_c': 17,
    'M_THREO__45__DS__45__ISO__45__CITRATE_c': 13,
    'M_THR_c': 17,
    'M_THZ__45__P_c': 4,
    'M_THZ_c': 4,
    'M_TRP_c': 17,
    'M_TYRAMINE_c': 1,
    'M_TYR_c': 17,
    'M_UDP__45__AA__45__GLUTAMATE_c': 3,
    'M_UDP__45__ACETYL__45__CARBOXYVINYL__45__GLUCOSAMINE_c': 3,
    'M_UDP__45__D__45__GALACTO__45__14__45__FURANOSE_c': 7,
    'M_UDP__45__MANNACA_c': 1,
    'M_UDP__45__MANNAC_c': 3,
    'M_UDP__45__N__45__ACETYLMURAMATE_c': 3,
    'M_UDP__45__N__45__ACETYL__45__D__45__GLUCOSAMINE_c': 3,
    'M_UDP_c': 10,
    'M_UMP_c': 10,
    'M_URACIL_c': 3,
    'M_URATE_c': 17,
    'M_UREA_c': 17,
    'M_URIDINE_c': 10,
    'M_UROCANATE_c': 4,
    'M_UTP_c': 10,
    'M_VAL_c': 17,
    'M_VITAMIN_D3_c': 17,
    'M_WATER_c': 17,
    'M_XANTHINE_c': 10,
    'M_XANTHOSINE__45__5__45__PHOSPHATE_c': 10,
    'M_XANTHOSINE_c': 10,
    'M_XTP_c': 7,
    'M_XYLITOL_c': 17,
    'M_XYLULOSE__45__5__45__PHOSPHATE_c': 16,
    'M_ZN__43__2_c': 17
 }


def test_m2m_cscope_call():
    """
    Test m2m addedvalue when called in terminal.
    """
    # RUN THE COMMAND
    inppath = 'metabolic_data/'
    respath = 'addedvalue_output/'
    if not os.path.exists(respath):
        os.makedirs(respath)
    with tarfile.open(inppath + 'toy_bact.tar.gz') as tar:
        tar.extractall(path=respath)
    subprocess.call([
        'm2m', 'iscope', '-n', respath + '/toy_bact', '-o',
        respath, '-s', inppath + '/seeds_toy.sbml',
        '-q'
    ])
    iscope_file = respath + 'indiv_scopes/rev_iscope.tsv'

    compounds = {}
    compounds_name = {}
    # ISCOPE ANALYSIS
    with open(iscope_file, 'r') as iscope_matrix:
        iscope_reader = csv.reader(iscope_matrix, delimiter='\t')
        for index, compound in enumerate(next(iscope_reader)):
            compounds_name[index] = compound
        for line in iscope_reader:
            for index in range(len(line)):
                if index not in compounds:
                    if line[index] == '1':
                        compounds[index] = 1
                    else:
                        compounds[index] = 0
                else:
                    if line[index] == '1':
                        compounds[index] +=1

    named_compounds = {compounds_name[index] : compounds[index] for index in compounds}

    for compound in EXPECTED_PRODUCED_COMPOUNDS:
        assert EXPECTED_PRODUCED_COMPOUNDS[compound] ==  named_compounds[compound]

    # clean
    shutil.rmtree(respath)

if __name__ == "__main__":
    test_m2m_cscope_call()