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
import networkx as nx
import os
import shutil
import subprocess
import sys
import time
import zipfile

from bubbletools import convert
from metage2metabo import utils
from metage2metabo.m2m_analysis.taxonomy import extract_taxa, get_taxon

# Deactivate clingo module to avoid issue like this one:
# https://github.com/Aluriak/PowerGrASP/issues/1
import clyngor
clyngor.deactivate_clingo_module()

logger = logging.getLogger(__name__)


def powergraph_analysis(gml_input_file_folder, output_folder, oog_jar=None, taxon_file=None, taxonomy_level="phylum"):
    """Run the graph compression and picture creation

    Args:
        gml_input_file_folder (str): path to the gml folder or the gml file
        output_folder (str): path to the output folder
        oog_jar (str): path to OOG jar file
        taxon_file (str): mpwt taxon file for species
        taxonomy_level (str): taxonomy level, must be: phylum, class, order, family, genus or species.
    """
    starttime = time.time()

    gml_paths = utils.file_or_folder(gml_input_file_folder)

    bbl_path = os.path.join(output_folder, 'bbl')
    svg_path = os.path.join(output_folder, 'svg')
    html_output = os.path.join(output_folder, 'html')

    if not utils.is_valid_dir(bbl_path):
        logger.critical("Impossible to access/create output directory " + bbl_path)
        sys.exit(1)

    if oog_jar:
        if not utils.is_valid_dir(svg_path):
            logger.critical("Impossible to access/create output directory " +  svg_path)
            sys.exit(1)

    if not utils.is_valid_dir(html_output):
        logger.critical("Impossible to access/create output directory " + html_output)
        sys.exit(1)

    # 26 colours from Alphabet project (minus white):
    # https://en.wikipedia.org/wiki/Help:Distinguishable_colors
    alphabet_project_distinct_hex_colors = ["#F0A3FF",
        "#0075DC", "#993F00", "#4C005C", "#191919",
        "#005C31", "#2BCE48", "#FFCC99", "#808080",
        "#94FFB5", "#8F7C00", "#9DCC00", "#C20088",
        "#003380", "#FFA405", "#FFA8BB", "#426600",
        "#FF0010", "#5EF1F2", "#00998F", "#E0FF66",
        "#740AFF", "#990000", "#FFFF80", "#FFFF00",
        "#FF5005"]

    # 269 colors from:
    # https://graphicdesign.stackexchange.com/a/3815
    hex_colors = ["#000000","#FFFF00","#1CE6FF","#FF34FF",
    "#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900",
    "#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87",
    "#5A0007","#809693","#FEFFE6","#1B4400","#4FC601","#3B5DFF",
    "#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0",
    "#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035",
    "#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F",
    "#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2",
    "#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF",
    "#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F",
    "#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C",
    "#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9",
    "#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF",
    "#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500",
    "#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700",
    "#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0",
    "#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3",
    "#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E",
    "#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800",
    "#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4",
    "#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3",
    "#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01",
    "#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D",
    "#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F",
    "#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405",
    "#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78",
    "#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4",
    "#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF",
    "#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702",
    "#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E",
    "#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0",
    "#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE",
    "#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2",
    "#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A",
    "#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183",
    "#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6",
    "#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03",
    "#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD",
    "#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E",
    "#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8",
    "#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188",
    "#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71",
    "#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66",
    "#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46",
    "#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C",
    "#92896B"]

    if taxon_file:
        taxonomy_output_file = os.path.join(output_folder, 'taxonomy_species.tsv')
        tree_output_file = os.path.join(output_folder, 'taxon_tree.txt')
        if not os.path.exists(taxonomy_output_file):
            extract_taxa(taxon_file, taxonomy_output_file, tree_output_file, taxonomy_level)

        taxon_species, all_taxons = get_taxon(taxonomy_output_file)

    for target_name in gml_paths:
        bbl_output = os.path.join(bbl_path, target_name + '.bbl')
        svg_file = os.path.join(svg_path, target_name + '.bbl.svg')

        html_target = os.path.join(html_output, target_name)
        if not utils.is_valid_dir(html_target):
            logger.critical("Impossible to access/create output directory " + html_target)
            sys.exit(1)

        # Compress gml file into a bbl file with PowerGrASP.
        gml_input_path = gml_paths[target_name]
        logger.info('######### Graph compression: ' + target_name + ' #########')
        compression(gml_input_path, bbl_output)
        logger.info('######### PowerGraph visualization: ' + target_name + ' #########')

        # Read gml file with networkx and extract the essential and alternative symbionts using the note of each node (organism).
        graph = nx.read_gml(gml_input_path)
        essentials = [organism for organism in graph.nodes if graph.nodes[organism]['note'] == 'ES']
        alternatives = [organism for organism in graph.nodes if graph.nodes[organism]['note'] == 'AS']
        if taxon_file:
            key_species = essentials + alternatives
            taxon_key_species = set([organism.split('__')[0] for organism in key_species])
            if len(set(all_taxons).intersection(taxon_key_species)) == 0:
                logger.critical('Difference of taxonomy level between gml file ('+gml_input_path+') compared to '+taxonomy_output_file+'.')
                sys.exit(1)

            # For each taxon in the key species, create a color.
            if len(taxon_key_species) <= len(alphabet_project_distinct_hex_colors):
                used_colors = alphabet_project_distinct_hex_colors
            elif len(hex_colors) >= len(taxon_key_species) > len(alphabet_project_distinct_hex_colors):
                used_colors = hex_colors
            else:
                logger.critical('Too many taxa in key species (superior to 269) and not enough colors in palette to color the powergraph.')
                sys.exit(1)

            taxon_colors = {}
            for index, taxon in enumerate(taxon_key_species):
                taxon_colors[taxon] = used_colors[index]

        bbl_to_html(bbl_output, html_target)
        if taxon_file:
            if os.path.exists(html_target +'_taxon'):
                shutil.rmtree(html_target +'_taxon')
            shutil.copytree(html_target, html_target +'_taxon')
            update_js_taxonomy(html_target +'_taxon', taxon_colors, essentials, alternatives)
            output_html_merged = os.path.join(html_output, target_name + '_powergraph_taxon.html')
            merge_html_css_js(html_target +'_taxon', output_html_merged)

        update_js(html_target, essentials, alternatives)
        output_html_merged = os.path.join(html_output, target_name + '_powergraph.html')
        merge_html_css_js(html_target, output_html_merged)

        if oog_jar:
            svg_file = os.path.join(svg_path, target_name + '.bbl.svg')
            if os.path.exists(svg_file):
                os.remove(svg_file)
            bbl_to_svg(oog_jar, bbl_output, svg_path)
            if taxon_file:
                taxonomy_svg_file = os.path.join(svg_path, target_name + '_taxon.bbl.svg')
                if os.path.exists(taxonomy_svg_file):
                    os.remove(taxonomy_svg_file)
                shutil.copyfile(svg_file, taxonomy_svg_file)
                update_svg_taxonomy(taxonomy_svg_file, taxon_colors)
            update_svg(svg_file, essentials, alternatives)

    logger.info(
        "--- Powergraph runtime %.2f seconds ---\n" % (time.time() - starttime))


def compression(gml_input, bbl_output):
    """Solution graph compression

    Args:
        gml_input (str): gml file
        bbl_output (str): bbl output file
    """
    starttime = time.time()
    with open('powergrasp.cfg', 'w') as config_file:
        config_file.write('[powergrasp options]\n')
        config_file.write('SHOW_STORY = no\n')

    import powergrasp
    from bubbletools import BubbleTree

    with open(bbl_output, 'w') as fd:
        for line in powergrasp.compress_by_cc(gml_input):
            fd.write(line + '\n')

    tree = BubbleTree.from_bubble_file(bbl_output)
    logger.info('Number of powernodes: ' + str(len([powernode for powernode in tree.powernodes()])))
    logger.info('Number of poweredges: ' + str(tree.edge_number()))

    os.remove('powergrasp.cfg')

    logger.info(
        'Compression runtime %.2f seconds ---\n' % (time.time() - starttime))


def check_oog_jar_file(oog_jar):
    """Check Oog jar file

    Args:
        oog_jar (str): path to oog jar file
    """
    if not os.path.isfile(oog_jar):
        sys.exit('Check Oog.jar: ' + oog_jar + ' is not an available file.')

    try:
        jarfile = zipfile.ZipFile(oog_jar, "r")
    except zipfile.BadZipFile:
        sys.exit('Check Oog.jar: ' + oog_jar + ' is not a valid .jar file (as it is not a correct zip file).')

    oog_class = None
    manifest_jar = None

    for filename in jarfile.namelist():
        if filename.endswith('Oog.class'):
            oog_class = True
        if filename.endswith('MANIFEST.MF'):
            manifest_jar = True

    jarfile.close()

    if oog_class and manifest_jar:
        return True
    elif manifest_jar:
        logger.info('Check Oog.jar: no correct Oog.class in jar file ' + oog_jar)
        return True
    else:
        sys.exit('Check Oog.jar: not a correct jar file ' + oog_jar)


def bbl_to_html(bbl_input, html_output):
    """Powergraph website creation.
    This create a folder with html/CSS/JS files. By using the index.html file in a browser, user can see the powergraph.

    Args:
        bbl_input (str): bbl input file
        html_output (str): html output file
    """
    logger.info('######### Creation of the powergraph website accessible at ' + html_output + ' #########')
    convert.bubble_to_js(bbl_input, html_output, width_as_cover=False)


def bbl_to_svg(oog_jar, bbl_input, svg_output):
    """Powergraph picture creation

    Args:
        oog_jar (str): path to oog jar file
        bbl_input (str): bbl input file
        svg_output (str): svg output file
    """
    if not shutil.which('java'):
        logger.critical('java is not in the Path, m2m_analysis option --oog can not work without it.')
        sys.exit(1)

    check_oog = check_oog_jar_file(oog_jar)

    if check_oog:
        logger.info('######### Creation of the powergraph svg accessible at ' + svg_output + ' #########')
        oog_cmds = ["java", "-jar", oog_jar, "-inputfiles=" + bbl_input, "-img", "-outputdir=" + svg_output]
        subproc = subprocess.Popen(oog_cmds)
        subproc.wait()


def update_js(html_output, essentials, alternatives):
    """Update graph.js to add colors for essential and alternative symbionts.

    Args:
        html_output (str): path to html folder (containing js subfolder with gaph.js)
        essentials (list): list of essential symbionts
        alternatives (list): list of alternative symbionts
    """
    # Add selector to color node according to the type if node:
    # circle #D41159 for essential symbiont
    # round rectangle #1A85FF for alternative symbiont
    selector_color = '''
    {
        selector: 'node[type="essential"]',
        css: {
            'background-color': '#D41159',
            'height': '30px',
        }
    },
    {
        selector: 'node[type="alternative"]',
        css: {
            'background-color': '#1A85FF',
            'width': '30px',
            'shape':'rectangle'
        }
    },
    '''

    js_folder = os.path.join(html_output, 'js')
    graph_js = os.path.join(js_folder, 'graph.js')
    new_graph_sj = ''
    with open(graph_js, 'r') as input_js:
        for line in input_js:
            if "data: { 'id'" in line:
                species_id = line.split("'id':")[1].split(',')[0].split("'")[1]
                if species_id in essentials:
                    line = line.replace(" } },", ", 'type': 'essential' } },")
                if species_id in alternatives:
                    line = line.replace(" } },", ", 'type': 'alternative' } },")
            new_graph_sj += line
            if 'style: [' in line:
                new_graph_sj += selector_color

    with open(graph_js, 'w') as input_js:
        input_js.write(new_graph_sj)


def update_js_taxonomy(html_output, taxon_colors, essentials, alternatives):
    """Update graph.js to add colors according to taxon.

    Args:
        html_output (str): path to html folder (containing js subfolder with gaph.js)
        taxon_colors (dict): dictionary {taxon_name: associated_color}
        essentials (list): list of essential symbionts
        alternatives (list): list of alternative symbionts
    """
    # Add selector to color node according to the type of node:
    # link a color to each taxon

    selector_color = ''
    for taxon in taxon_colors:
        taxon_color = taxon_colors[taxon]
        selector_color += '''
        {
            selector: 'node[type="'''+taxon+'''"]',
            css: {
                'background-color': "'''+taxon_color+'''",
            }
        },
        '''

    js_folder = os.path.join(html_output, 'js')
    graph_js = os.path.join(js_folder, 'graph.js')
    new_graph_sj = ''
    with open(graph_js, 'r') as input_js:
        for line in input_js:
            if "data: { 'id'" in line:
                species_names = line.split("'id':")[1].split(',')[0].strip("'| ")
                species_taxon_id = species_names.split('__')[0]
                if species_taxon_id in taxon_colors:
                    if species_names in essentials:
                        line = line.replace(" } },", ", 'type': '"+species_taxon_id+"' } },")
                    if species_names in alternatives:
                        line = line.replace(" } },", ", 'type': '"+species_taxon_id+"' }, css: {'shape': 'rectangle'} },")

            new_graph_sj += line
            if 'style: [' in line:
                new_graph_sj += selector_color

    with open(graph_js, 'w') as input_js:
        input_js.write(new_graph_sj)


def update_svg(svg_file, essentials, alternatives):
    """Update svg file to add colors for essential and alternative symbionts.

    Args:
        svg_file (str): path to svg file
        essentials (list): list of essential symbionts
        alternatives (list): list of alternative symbionts
    """
    new_svg = []
    previous_line = ''
    with open(svg_file, 'r') as input_js:
        for line in input_js:
            if '<text' in line:
                species_id = line.split('>')[1].split('<')[0]
                if species_id in essentials:
                    new_color = ' style="stroke:none;fill:#D41159"'
                    line_before_style = line.split(' style="')[0]
                    line_after_style = '"'.join(line.split(' style="')[1].split('"')[1:])
                    line = line_before_style + new_color + line_after_style
                    previous_line_before_style = previous_line.split(' style="')[0]
                    previous_line_after_style = '"'.join(previous_line.split(' style="')[1].split('"')[1:])
                    previous_line = previous_line_before_style + new_color + previous_line_after_style
                    previous_line = previous_line.replace('circle', 'ellipse rx="5" ry="5"')
                    new_svg[-1] = previous_line
                if species_id in alternatives:
                    new_color = ' style="stroke:none;fill:#1A85FF"'
                    line_before_style = line.split(' style="')[0]
                    line_after_style = '"'.join(line.split(' style="')[1].split('"')[1:])
                    line = line_before_style + new_color + line_after_style
                    previous_line_before_style = previous_line.split(' style="')[0]
                    previous_line_after_style = '"'.join(previous_line.split(' style="')[1].split('"')[1:])
                    previous_line = previous_line_before_style + new_color + previous_line_after_style
                    previous_line = previous_line.replace('circle', 'ellipse rx="5" ry="5"')
                    new_svg[-1] = previous_line
            new_svg.append(line)
            previous_line = line

    with open(svg_file, 'w') as input_js:
        input_js.write(''.join(new_svg))


def update_svg_taxonomy(svg_file, taxon_colors):
    """Update svg file to add colors for each taxon

    Args:
        svg_file (str): path to svg file
        taxon_colors (dict): dictionary {taxon_name: associated_color}
    """
    new_svg = []
    previous_line = ''
    with open(svg_file, 'r') as input_js:
        for line in input_js:
            if '<text' in line:
                species_taxon_id = line.split('>')[1].split('<')[0].split('__')[0]
                if species_taxon_id in taxon_colors:
                    taxon_color = taxon_colors[species_taxon_id]
                    new_color = ' style="stroke:none;fill:{0};"'.format(taxon_color)
                    line_before_style = line.split(' style="')[0]
                    line_after_style = '"'.join(line.split(' style="')[1].split('"')[1:])
                    line = line_before_style + new_color + line_after_style
                    previous_line_before = previous_line.split(' style="')[0]
                    previous_line_after_style = '"'.join(previous_line.split(' style="')[1].split('"')[1:])
                    previous_line = previous_line_before + new_color + previous_line_after_style
                    new_svg[-1] = previous_line
            new_svg.append(line)
            previous_line = line

    with open(svg_file, 'w') as input_js:
        input_js.write(''.join(new_svg))


def merge_html_css_js(html_output, merged_html_path):
    """Merge HTML/CSS/JS files into one HTML file

    Args:
        html_output (str): path to html folder (containing css, html and css files)
        merged_html_path (str): path to the output merged html file
    """
    index_html = os.path.join(html_output, 'index.html')
    style_css = os.path.join(html_output, 'style.css')
    js_folder = os.path.join(html_output, 'js')
    graph_js = os.path.join(js_folder, 'graph.js')
    cytoscape_min_js = os.path.join(js_folder, 'cytoscape.min.js')
    cytoscape_cose_bilkent_js = os.path.join(js_folder, 'cytoscape-cose-bilkent.js')

    output_html = os.path.join(merged_html_path)

    with open(style_css, 'r') as input_css:
        css_str =  input_css.read()

    with open(graph_js, 'r') as input_graph_js:
        graph_js_str =  input_graph_js.read()

    with open(cytoscape_min_js, 'r') as input_cytoscape_min_js:
        cytoscape_min_js_str =  input_cytoscape_min_js.read()

    with open(cytoscape_cose_bilkent_js, 'r') as input_cytoscape_cose_bilkent_js:
        cytoscape_cose_bilkent_js_str =  input_cytoscape_cose_bilkent_js.read()

    new_html_str = ''
    line_before = ''
    with open(index_html, 'r') as input_html_file:
        for line in input_html_file:
            if '<head>' in line_before:
                new_html_str += '        <style>'
                new_html_str += css_str
                new_html_str += '        </style>'
            if 'meta name="viewport"' in line_before:
                new_html_str += '    <script>'
                new_html_str += cytoscape_min_js_str
                new_html_str += cytoscape_cose_bilkent_js_str
                new_html_str += '    </script>'
            if 'div id="cy"' in line_before:
                new_html_str += '    <script>'
                new_html_str += graph_js_str
                new_html_str += '    </script>'
            if 'script' not in line and 'link' not in line:
                new_html_str +=  line
            line_before = line

    with open(output_html, 'w') as output_hml_file:
        output_hml_file.write(new_html_str)
