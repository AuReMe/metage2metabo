============
m2m_analysis
============
m2m_analysis is a supplementary command installed by metage2metabo. It is separated from m2m because it has heavy dependencies. Also, m2m_analysis steps can be time consuming.

Requirements
------------

m2m_analysis needs:

* `Oog Power Graph Command line tool <http://www.biotec.tu-dresden.de/research/schroeder/powergraphs/download-command-line-tool.html>`__: used to create a svg file of the powergraph. It is a jar file (compiled for Java 6), so you need at least Java 6.

* Some python packages:

    * `networkx <https://github.com/networkx/networkx>`__: to create graph from miscoto results
    * `ete3 <https://github.com/etetoolkit/ete>`__: to add taxonomy information on the graph if you used mpwt taxon file
    * `powergrasp <https://github.com/Aluriak/PowerGrASP>`__: to compress networkx graph (which required graphviz)

Requirements
------------

m2m_analysis goes deeper in the analysis compare to m2m. In m2m_analysis, the enumeration of all solution is computed, this step is far more time consuming than the others (union or intersection).

This first step is done by ``m2m_analysis enum``. It will enumerate all the possible solutions, select the optimal ones and then create the json file containing all the optimal solutions.

In a second step (``m2m_analysis graph``), the optimal solutions from the enumeration are stored in a graph. The nodes of this graph are each organism present in at least one solution. An edge connects two nodes, only if the two organisms (represented by the node) co-occur in at least one of the enumerated communities.

The last step (``m2m_analysis powergraph``) compresses the graph into a power graph (in bbl format). Then it creates a svg picture of this power graph.

m2m_analysis enum
-----------------
``m2m_analysis enum`` runs miscoto with enumeration on a set of target.

It uses the following mandatory inputs (run ``m2m_analysis enum --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-t file/directory      targets SBML file or folder containing multiple targets SBML files
-o directory           output directory for results
-m file                host metabolic network in SBML

m2m_analysis graph
------------------
``m2m_analysis graph`` creates the graph containing the solutions.

It uses the following mandatory inputs (run ``m2m_analysis graph --help`` for optional arguments):

-j file/directory      directory of miscoto output JSONs
-t file/directory      targets SBML file or folder containing multiple targets SBML files
-o directory           output directory for results
--taxon file           mpwt taxon file

m2m_analysis powergraph
-----------------------
``m2m_analysis powergraph`` compresses the graph and create a svg picture.

It uses the following mandatory inputs (run ``m2m_analysis powergraph --help`` for optional arguments):

--oog file             Oog jar file
-g file/directory      directory of GML files or a GML file
-o directory           output directory for results

m2m_analysis workflow
---------------------
``m2m_analysis workflow`` runs the all m2m_analysis workflow.

It uses the following mandatory inputs (run ``m2m_analysis workflow --help`` for optional arguments):

-n directory           directory of metabolic networks, 
                        in SBML format
-s file                seeds SBML file
-t file/directory      targets SBML file or folder containing multiple targets SBML files
-o directory           output directory for results
-m file                host metabolic network in SBML
--oog file             Oog jar file
--taxon file           mpwt taxon file