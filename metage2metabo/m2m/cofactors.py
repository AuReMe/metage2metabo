# Copyright (C) 2019-2021 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import csv
import logging
import os
import sys
import tempfile
import time

from metage2metabo import utils

from moped import Model
from padmet.utils.sbmlPlugin import convert_to_coded_id
from multiprocessing import Pool


logger = logging.getLogger(__name__)
logging.getLogger("moped").setLevel(logging.ERROR)
logging.getLogger("cobra").setLevel(logging.CRITICAL)


def create_mock_cofactors(input_sbml_path, pair_cofactors, output_sbml_path):
    sbml_basename = utils.get_basename(input_sbml_path)

    m = Model()
    m.read_from_sbml(input_sbml_path)

    # Add pair of cofactors.
    m.cofactor_pairs = pair_cofactors

    # Use moped duplication of cofactors.
    m.cofactor_duplication()

    # Extract reactions and metabolites added by moped duplication.
    rxn_cofactors = {}
    metabolite_cofactors = []
    reversibility_reactions = {}
    for reaction in m.reactions:
        products = []
        reactants = []
        if '__cof__' in reaction:
            for metabolite in m.reactions[reaction].stoichiometries:
                if m.reactions[reaction].stoichiometries[metabolite] > 0:
                    products.append(convert_to_coded_id(metabolite))
                elif m.reactions[reaction].stoichiometries[metabolite] < 0:
                    reactants.append(convert_to_coded_id(metabolite))
                if '__cof__' in metabolite:
                    metabolite_cofactors.append(convert_to_coded_id(metabolite))
            rxn_cofactors[reaction] = (products, reactants)
            reversibility_reactions[reaction] = m.reactions[reaction].reversible

    metabolite_cofactors = set(metabolite_cofactors)

    # Add them to the sbml in a very ugly way.
    # If I believed in hell, I think that lines of code would guarantee me an access to that place.
    # This add the cofactor metabolites at the end of the listOfSpecies and the cofactors reactions at the end of the listOfReactions.
    sbml_text = ''

    metabolite_sbml = '			<species id="{0}"/>'
    reaction_sbml = '      <reaction id="{0}"  reversible="{1}" >'
    species_ref_sbml = '          <speciesReference species="{0}" stoichiometry="1" constant="true"/>'

    with open(input_sbml_path, 'r') as sbml_input_file:
        for line in sbml_input_file:
            if '/listOfSpecies' in line:
                for metabolite_cof in metabolite_cofactors:
                    sbml_text += metabolite_sbml.format(metabolite_cof)+'\n'
            if '/listOfReactions' in line:
                for reaction in rxn_cofactors:
                    sbml_text += reaction_sbml.format(convert_to_coded_id(reaction), str(reversibility_reactions[reaction]).lower())+'\n'
                    sbml_text += '        <listOfReactants>'+'\n'
                    for reactant in rxn_cofactors[reaction][0]:
                        sbml_text += species_ref_sbml.format(reactant)+'\n'
                    sbml_text += '        </listOfReactants>'+'\n'
                    sbml_text += '        <listOfProducts>'+'\n'
                    for product in rxn_cofactors[reaction][1]:
                        sbml_text += species_ref_sbml.format(product)+'\n'
                    sbml_text += '        </listOfProducts>\n      </reaction>\n'
            sbml_text += line

    # Write the modified sbml to the output path.
    with open(output_sbml_path, 'w') as sbml_output_file:
        sbml_output_file.write(sbml_text)

    logger.info('Add {0} mocked metabolites and {1} mocked reactions for {2}.'.format(len(metabolite_cofactors), len(rxn_cofactors), sbml_basename))


def cofactors(sbml_dir, cofactors_file, output_folder, cpu_number=1):
    cofactors_output = os.path.join(output_folder, 'mocked_sbmls')

    if not utils.is_valid_dir(cofactors_output):
        logger.critical('Impossible to access/create output directory')
        sys.exit(1)

    cofactor_pool = Pool(cpu_number)

    all_files = [
        f for f in os.listdir(sbml_dir)
        if os.path.isfile(os.path.join(sbml_dir, f)) and utils.get_extension(
            os.path.join(sbml_dir, f)).lower() in ['xml', 'sbml']
    ]

    pair_cofactors = {}
    with open(cofactors_file, 'r') as cofactor_open_file:
        csvreader = csv.reader(cofactor_open_file, delimiter='\t')
        # Avoid header.
        next(csvreader)
        for line in csvreader:
            pair_cofactors[line[0]] = line[1]

    multiprocessing_cofactors = []
    for f in all_files:
        sbml_path = os.path.join(sbml_dir, f)
        output_file_path = os.path.join(cofactors_output, f)
        multiprocessing_cofactors.append((sbml_path, pair_cofactors, output_file_path))

    cofactor_pool.starmap(create_mock_cofactors, multiprocessing_cofactors)

    cofactor_pool.close()
    cofactor_pool.join()
