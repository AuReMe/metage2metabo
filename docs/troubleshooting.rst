===============
Troubleshooting
===============

Several potential fix for issues we encounter.

`m2m_analysis`: "Can't connect to X11 window server using..."
-------------------------------------------------------------

Issue encountered when using `m2m_analysis` on a HPC:

::
    
    <II> Create SVG ... Exception in thread "main" java.awt.AWTError: Can't connect to X11 window server using ':0' as the value of the DISPLAY variable.
        at sun.awt.X11GraphicsEnvironment.initDisplay(Native Method)
        at sun.awt.X11GraphicsEnvironment.access$200(X11GraphicsEnvironment.java:65)
        at sun.awt.X11GraphicsEnvironment$1.run(X11GraphicsEnvironment.java:115)
        at java.security.AccessController.doPrivileged(Native Method)
        at sun.awt.X11GraphicsEnvironment.<clinit>(X11GraphicsEnvironment.java:74)
        at java.lang.Class.forName0(Native Method)
        at java.lang.Class.forName(Class.java:264)
        at java.awt.GraphicsEnvironment.createGE(GraphicsEnvironment.java:103)
        at java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment(GraphicsEnvironment.java:82)
        at java.awt.image.BufferedImage.createGraphics(BufferedImage.java:1181)
        at org.apache.batik.svggen.SVGGraphics2D.<init>(Unknown Source)
        at org.apache.batik.svggen.SVGGraphics2D.<init>(Unknown Source)
        at org.mattlab.eaglevista.ui.GraphImageFactory.getSVGGraphicsOfGraph(GraphImageFactory.java:157)
        at org.mattlab.eaglevista.ui.GraphImageFactory.createImage(GraphImageFactory.java:80)
        at org.mattlab.eaglevista.ui.EvCommandLineInterpreter.runImpl(EvCommandLineInterpreter.java:466)
        at org.mattlab.eaglevista.ui.EvCommandLineInterpreter.main(EvCommandLineInterpreter.java:81)


One potential fix is to unset DISPLAY before launching m2m:

::

    unset DISPLAY
    m2m_analysys workflow...


UnicodeDecodeError when running m2m_analysis with Singularity
-------------------------------------------------------------

Error encountered when using `m2m_analysis` in a Singularity:

::

    ######### Creation of the powergraph website accessible at results_m2m_analysis_scfa/html/metabolic_targets_scfa #########
    Traceback (most recent call last):
    File "/usr/local/bin/m2m_analysis", line 8, in <module>
        sys.exit(main())
    File "/usr/local/lib/python3.6/dist-packages/metage2metabo/__main_analysis__.py", line 289, in main
        args.oog, new_arg_modelhost, args.level)
    File "/usr/local/lib/python3.6/dist-packages/metage2metabo/__main_analysis__.py", line 303, in main_analysis_workflow
        run_analysis_workflow(*allargs)
    File "/usr/local/lib/python3.6/dist-packages/metage2metabo/m2m_analysis/m2m_analysis_workflow.py", line 48, in run_analysis_workflow
        powergraph_analysis(gml_output, output_dir, oog_jar, taxon_file, taxonomy_level)
    File "/usr/local/lib/python3.6/dist-packages/metage2metabo/m2m_analysis/graph_compression.py", line 174, in powergraph_analysis
        merge_html_css_js(html_target +'_taxon', output_html_merged)
    File "/usr/local/lib/python3.6/dist-packages/metage2metabo/m2m_analysis/graph_compression.py", line 473, in merge_html_css_js
        cytoscape_min_js_str =  input_cytoscape_min_js.read()
    File "/usr/lib/python3.6/encodings/ascii.py", line 26, in decode
        return codecs.ascii_decode(input, self.errors)[0]
    UnicodeDecodeError: 'ascii' codec can't decode byte 0xe2 in position 234868: ordinal not in range(128)

This is an issue with the encodign standard used. A solution is to set the local codec to `LC_ALL=C.UTF-8`:

::

    LC_ALL=C.UTF-8 && srun singularity exec -B....

