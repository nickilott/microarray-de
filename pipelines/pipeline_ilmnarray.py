################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
Pipeline pipeline_ilmnarray
===========================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline for the analysis of illumina gene expression.

Overview
========

There are some general goals when conducting an expression profiling experiment. Which 
genes are differentially expressed between conditions and which pathways are those
genes associated with. Background correction, normalisation and differential expression
analysis are standard steps in the processing of these data and are therefore amenable to
automation. The pipeline performs these steps along with pathways analysis for the rapid 
and automated analysis and reporting of illumina beadarrays.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

The input is data that has been processed by BeadStudio software. This is in part because
most of the arrays that we receive at the NDM are processed in this way as they are run
at the Wellcome Trust Centre for Human Genetics (WTCHG).

The BeadStudio output must be named sample_probe_profile.txt.

The headers in the file must conform to the following:

    <tissue>-<treatment>.<replicate>

for example:

    treg-il10neg.R1

This is to allow the pipeline to pick up the different conditions that are to be tested against each other.


In order for the pipeline to run pathways enrichment analysis there needs to be an additional file in 
the working directory. This file maps gene names (present in the SampleProbeProfile.txt file) to pathway
identifiers and descriptions. It does not matter which database is used but the gene names MUST match.

The input file looks like:


+--------+-------+--------+-----------------------------------------------+--------+
|ontology|gene_id|kegg_ID |kegg_name                                      |evidence|
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |A2M    |hsa04610|"Complement and coagulation cascades"          |        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |A4GALT |hsa00603|"Glycosphingolipid biosynthesis - globo series"|        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |A4GALT |hsa01100|"Metabolic pathways"                           |        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |AAAS   |hsa03013|"RNA transport"                                |        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |AACS   |hsa00650|"Butanoate metabolism"                         |        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |AADAT  |hsa00300|"Lysine biosynthesis"                          |        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |AADAT  |hsa00310|"Lysine degradation"                           |        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |AADAT  |hsa00380|"Tryptophan metabolism"                        |        |
+--------+-------+--------+-----------------------------------------------+--------+
|kegg    |AADAT  |hsa01100|"Metabolic pathways"                           |        |
+--------+-------+--------+-----------------------------------------------+--------+
 
where the gene ids are gene names and the pathwayds are from the KEGG pathways database.

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following R packages to be installed:


+--------------------+------------------------------------------------+
|*Program*           |*Purpose*                                       |
+--------------------+------------------------------------------------+
| LIMMA              | Normalization and differential expression      |
+--------------------+------------------------------------------------+
| pheatmap           | Drawing pretty heatmaps                        |
+--------------------+------------------------------------------------+
| pvclust            | Clustering using bootstrap resampling          |
+--------------------+------------------------------------------------+
| ggplot2            | plotting library                               |
+--------------------+------------------------------------------------+
| affy               | affy for normalisation                         |
+--------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

This contains the results of differential expression analysis and pathway enrichment analysis. In addition
to this, automated reports are generated.

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3
from rpy2.robjects import r as R

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import PipelineIlmnArray

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters( 
    [ "pipeline.ini" ] )

PARAMS = P.PARAMS

#PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
 #                                      "pipeline_annotations.py" )

###################################################################
###################################################################
###################################################################
def connect():
    '''
    connect to database.
    Use this method to connect to additional databases.
    Returns a database connection.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    # statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    # cc = dbh.cursor()
    # cc.execute( statement )
    # cc.close()

    return dbh

###################################################################
###################################################################
###################################################################
## worker tasks
###################################################################
CONTROL = {1: "control_probe_profile.txt",
           0: None}

@follows(mkdir("qc.dir"))
@split( ["sample_probe_profile.txt", CONTROL[PARAMS.get("norm_use_control")]],
        ["sample_probe_profile.matrix",
         "qc.dir/postprocess_density.pdf",
         "qc.dir/postprocess_boxplot.pdf",
         "qc.dir/mean_variance.pdf"])
def normaliseAndQc( infile, outfiles ):
    '''
    normalise the data and output qc plots.
    The main output is sample_probe_profile.matrix
    which is the normalised array data.
    '''
    # do not have this functionality (I think)
    to_cluster = False

    PipelineIlmnArray.normaliseAndQc(infile, 
                                     outfiles,
                                     PARAMS["general_rdir"],
                                     platform = PARAMS.get("platform"),
                                     filter_above = PARAMS["norm_filter_above"],
                                     filter_over = PARAMS["norm_filter_over"],
                                     preprocess_method = PARAMS["norm_preprocess_method"],
                                     preprocess = PARAMS["norm_preprocess"],
                                     method = PARAMS["norm_method"])

###################################################################
###################################################################
###################################################################
@follows(mkdir("probe_info.dir"))
@merge(["sample_probe_profile.txt", normaliseAndQc],
       "probe_info.dir/probe_info.tsv")
def buildProbeInfo(infiles, outfile):
    '''
    build a file with the total number of probes in etc
    '''
    to_cluster = False

    PipelineIlmnArray.buildProbeInfo(infiles,
                                     outfile)
        
###################################################################
###################################################################
###################################################################
@transform(buildProbeInfo, suffix(".tsv"), ".load")
def loadProbeInfo(infile, outfile):
    '''
    load probe_info file
    '''
    to_cluster = False
    
    tablename = P.toTable(outfile)
    statement = """/usr/local/bin/python2.7 %s/csv2db.py \
                   -t %s \
                   --log=%s.log \
                   < %s > %s""" % (PARAMS["general_cgatscriptsdir"], 
                                   tablename, 
                                   outfile, 
                                   infile, 
                                   outfile) 
    os.system(statement)

###################################################################
###################################################################
###################################################################
@transform("sample_probe_profile.txt", 
           suffix("sample_probe_profile.txt"),
           "probe2gene.map")
def buildProbe2GeneMap(infile, outfile):
    '''
    build file mapping probe ID to gene name
    '''
    to_cluster = False
    platform = PARAMS.get("platform")
    PipelineIlmnArray.buildProbe2GeneMap(infile, 
                                         outfile,
                                         PARAMS,
                                         platform = platform)

###################################################################
###################################################################
###################################################################
@transform(buildProbe2GeneMap, suffix(".map"), ".map.load")
def loadProbe2GeneMap(infile, outfile):
    '''
    load probe mapping to gene names
    '''
    to_cluster = False
    tablename = P.toTable(outfile)
    statement = """/usr/local/bin/python2.7 %s/csv2db.py \
                   -t %s \
                   --log=%s.log \
                   < %s > %s""" % (PARAMS["general_cgatscriptsdir"],
                                   tablename,
                                   outfile,
                                   infile,
                                   outfile)
                                   
    os.system(statement)

###################################################################
###################################################################
###################################################################
if normaliseAndQc:
    NORM_TARGET = "sample_probe_profile.matrix"

@follows(normaliseAndQc, mkdir("cluster.dir"))
@transform(NORM_TARGET, regex("(\S+).matrix"), r"cluster.dir/\1.pvclust.pdf")
def clusterSamplesOnExpressionValues(infile, outfile):
    '''
    use pvclust to cluster samples on expression values
    '''
    to_cluster = False
    PipelineIlmnArray.clusterSamplesOnExpressionValues(infile,
                                                       outfile,
                                                       dist_method = PARAMS["cluster_dist_method"],
                                                       linkage_method = PARAMS["cluster_linkage_method"],
                                                       nboot = int(PARAMS["cluster_nboot"]))

###################################################################
###################################################################
###################################################################
@follows(normaliseAndQc, mkdir("cluster.dir"))
@transform(NORM_TARGET, regex("(\S+).matrix"), r"cluster.dir/\1.pca.pdf")
def pcaSamplesOnExpressionValues(infile, outfile):
    '''
    use PCA analysis to cluster samples on expression values
    '''
    to_cluster = False
    PipelineIlmnArray.pcaSamplesOnExpressionValues(infile,
                                                   outfile)

###################################################################
###################################################################
###################################################################
@follows(normaliseAndQc, mkdir("cluster.dir"))
@transform(NORM_TARGET, regex("(\S+).matrix"), r"cluster.dir/\1.correlation.pdf")
def correlateSamplesOnExpressionValues(infile, outfile):
    '''
    correlate samples and draw a heatmap
    '''
    to_cluster = False
    PipelineIlmnArray.correlateSamplesOnExpressionValues(infile,
                                                   outfile)


###################################################################
###################################################################
###################################################################
## intermediate target
###################################################################
@follows( normaliseAndQc,
          buildProbeInfo,
          clusterSamplesOnExpressionValues,
          pcaSamplesOnExpressionValues,
          correlateSamplesOnExpressionValues)
def QC(): pass

###################################################################
###################################################################
###################################################################
# Differential expression analysis
###################################################################
###################################################################
###################################################################
@follows(normaliseAndQc, mkdir("differential_expression.dir"))
@transform(NORM_TARGET, regex("(\S+).matrix"), r"differential_expression.dir/\1.design")
def buildDesignFile(infile, outfile):
    '''
    use the conditions encoded in the header line to determine
    design of the experiment.
    '''
    to_cluster = False
    PipelineIlmnArray.buildDesignFile(infile, 
                                      outfile, 
                                      paired = PARAMS.get("paired"))

###################################################################
###################################################################
###################################################################
@follows(normaliseAndQc)
@split([NORM_TARGET, buildDesignFile, buildProbe2GeneMap], "differential_expression.dir/*.result")
def runLimma(infiles, outfiles):
    '''
    run limma for differential expression analysis
    '''
    to_cluster = False
    contrasts = PipelineIlmnArray.buildContrasts(infiles[1])
    contrasts = ",".join(contrasts)
    PipelineIlmnArray.runLimma(infiles, outfiles, 
                               contrasts, 
                               adjust_method = PARAMS["differential_expression_adjust_method"],
                               rdir = PARAMS["rdir"])

###################################################################
###################################################################
###################################################################
@transform(runLimma, suffix(".result"), ".result.load")
def loadLimma(infile, outfile):
    '''
    load limma results
    '''
    to_cluster = False
    tablename = P.toTable(outfile)
    statement = """/usr/local/bin/python2.7 %s/csv2db.py \
                   -t %s \
                   --log=%s.log \
                   < %s > %s""" % (PARAMS["general_cgatscriptsdir"],
                                   tablename,
                                   outfile,
                                   infile,
                                   outfile)
    os.system(statement)

###################################################################
###################################################################
###################################################################
@transform(runLimma, suffix(".result"), ".ma.pdf")    
def MAPlotResults(infile, outfile):
    '''
    MA plot the results
    '''
    to_cluster = False
    p_threshold = PARAMS["differential_expression_padj"]
    fc_threshold = PARAMS["differential_expression_l2fold"]
    PipelineIlmnArray.MAPlotResults(infile, outfile,
                                    rdir = PARAMS["rdir"]
                                    , p_threshold = p_threshold
                                    , fc_threshold = fc_threshold)
    
###################################################################
###################################################################
###################################################################
@transform(runLimma, suffix(".result"), ".volcano.pdf")    
def volcanoPlotResults(infile, outfile):
    '''
    MA plot the results
    '''
    to_cluster = False
    p_threshold = PARAMS["differential_expression_padj"]
    fc_threshold = PARAMS["differential_expression_l2fold"]
    PipelineIlmnArray.volcanoPlotResults(infile, outfile, 
                                         fc_threshold = fc_threshold,
                                         p_threshold = p_threshold,
                                         rdir = PARAMS["rdir"])
    
###################################################################
###################################################################
###################################################################
@merge([loadLimma, loadProbe2GeneMap], "differential_expression.dir/differentially_expressed.tsv")
def buildDifferentiallyExpressedGeneList(infiles, outfile):
    '''
    build an annotated set of probes on the array:
    0 = no change
    1 = up
    2 - down
    according to parameters specified in the .ini file
    '''
    to_cluster = False
    dbh = connect()
    tablenames = [P.toTable(inf) for inf in infiles]
    probe2gene_table = tablenames[-1]
    result_tablenames = tablenames[:-1]
    PipelineIlmnArray.buildDifferentiallyExpressedGeneList(infiles, 
                                                           outfile,
                                                           dbh,
                                                           result_tablenames,
                                                           probe2gene_table,
                                                           l2fold = PARAMS["differential_expression_l2fold"],
                                                           padj = PARAMS["differential_expression_padj"])

###################################################################
###################################################################
###################################################################
@transform(buildDifferentiallyExpressedGeneList, suffix(".tsv"), ".load")
def loadDifferentiallyExpressedGeneList(infile, outfile):
    '''
    load differentially expressed genes
    '''
    to_cluster = False
    tablename = P.toTable(outfile)
    statement = """/usr/local/bin/python2.7 %s/csv2db.py \
                   -t %s \
                   --log=%s.log \
                   < %s > %s""" % (PARAMS["general_cgatscriptsdir"],
                                   tablename,
                                   outfile,
                                   infile,
                                   outfile)

    os.system(statement)

###################################################################
###################################################################
###################################################################
@split(loadDifferentiallyExpressedGeneList, "differential_expression.dir/*.venn.pdf")
def vennOverlap(infile, outfiles):
    '''
    venn the overlap of gene lists from analysis
    '''
    dbh = connect()
    PipelineIlmnArray.vennOverlap(dbh, outfiles)



###################################################################
###################################################################
###################################################################
@transform(loadDifferentiallyExpressedGeneList, 
           suffix(".load"), 
           add_inputs(NORM_TARGET),
           ".heatmap.pdf")
def heatmapDifferentiallyExpressedGeneList(infiles, outfile):
    '''
    heatmap the expression values for differentially expressed
    genes
    '''
    to_cluster = False
    dbh = connect()
    diff_table = P.toTable(infiles[0])
    data_matrix = infiles[1]
    PipelineIlmnArray.heatmapDifferentiallyExpressedGeneList(dbh,
                                                             diff_table,
                                                             data_matrix,
                                                             outfile,
                                                             rdir = PARAMS["rdir"])

###################################################################
###################################################################
###################################################################
@follows(loadLimma,
         MAPlotResults,
         volcanoPlotResults,
         loadDifferentiallyExpressedGeneList,
         heatmapDifferentiallyExpressedGeneList,
         vennOverlap)
          

def differential_expression():
    pass

###################################################################
###################################################################
###################################################################
# Pathways analysis using KEGG annotations
###################################################################
###################################################################
###################################################################
@follows(mkdir("pathways.dir"))
@transform(buildDifferentiallyExpressedGeneList, regex("(\S+)/(\S+).tsv"), r"pathways.dir/\2.background")
def buildBackgroundFile(infile, outfile):
    '''
    build background set of genes - this is the list that
    was tested for differential expression
    '''
    to_cluster = False
    PipelineIlmnArray.buildBackgroundFile(infile, outfile)

###################################################################
###################################################################
###################################################################
@follows(mkdir("pathways.dir"))
@transform(buildDifferentiallyExpressedGeneList, regex("(\S+)/(\S+).tsv"), r"pathways.dir/\2.foreground")
def buildForegroundFile(infile, outfile):
    '''
    foreground file contains all comparisons and coded as:
    1 = differentially expressed (up or down)
    0 = not differentially expressed2
    '''
    to_cluster = False
    PipelineIlmnArray.buildForegroundFile(infile, outfile)

###################################################################
###################################################################
###################################################################
@split([buildBackgroundFile, buildForegroundFile, PARAMS["pathways_gene2pathway"]]
       , "pathways.dir/*.overall")
def runGO(infiles, outfiles):
    '''
    run go analysis 
    '''
    to_cluster = False
    background, genes, gene2pathway = infiles[0], infiles[1], infiles[2]

    statement = """/usr/local/bin/python2.7 %s/runGO.py \
                   --genes=%s \
                   --background=%s \
                   --filename-input=%s \
                   -q BH \
                   --fdr \
                   --output-filename-pattern="pathways.dir/%%(set)s.%%(go)s.%%(section)s" \
                   > pathways.dir/pathways.log  \
                """ % (PARAMS["general_cgatscriptsdir"], genes, background, gene2pathway)
    os.system(statement)

###################################################################
###################################################################
###################################################################
@jobs_limit(1, "R")
@transform(runGO, suffix(".overall"), ".barplot.pdf")
def barplotGO(infile, outfile):
    '''
    produce a barplot of significant results from pathways
    analysis
    '''
    to_cluster = False
    PipelineIlmnArray.barplotGO(infile, outfile)

###################################################################
###################################################################
###################################################################
@transform(runGO, regex("(\S+)/(\S+).overall"),
           add_inputs(PARAMS["pathways_gene2pathway"]),
           r"\1/\2.genes")
def buildPathwayGenes(infiles, outfile):
    '''
    build sets of genes that are differentially expressed and
    fall into the pathways identified in the top 10 lists
    '''
    goresult, gene2pathways = infiles[0], infiles[1]
    if "_up" in os.path.basename(goresult):
        contrast = os.path.basename(goresult).replace("_up", ",").split(",")[0]
    if "_down" in os.path.basename(goresult):
        contrast = os.path.basename(goresult).replace("_down", ",").split(",")[0]
    dbh = connect()
    result_table = contrast
    
    PipelineIlmnArray.buildPathwayGenes(goresult, 
                                        gene2pathways, 
                                        contrast, 
                                        dbh, 
                                        result_table, 
                                        outfile)

###################################################################
###################################################################
###################################################################
@transform(buildPathwayGenes, suffix(".genes"), ".genes.load")
def loadPathwayGenes(infile, outfile):
    '''
    load pathway genes
    '''
    tablename = P.toTable(outfile)
    statement = """/usr/local/bin/python2.7 %s/csv2db.py \
                   -t %s --ignore-empty \
                   < %s > %s \
                """ % (PARAMS["general_cgatscriptsdir"], tablename, infile, outfile)
    os.system(statement)

###################################################################
###################################################################
###################################################################
@jobs_limit(1, "R")
@transform(buildPathwayGenes, suffix(".genes"), ".genes.plots")
def plotPathwayGenes(infile, outfile):
    '''
    plot differentially expressed genes in identified 
    pathways
    creates sentinel file as many plots are produced - i.e
    one for each pathway
    '''
    PipelineIlmnArray.plotPathwayGenes(infile, outfile)


@follows(barplotGO
         , loadPathwayGenes
         )
def pathways():
    pass


## primary targets
###################################################################
@follows(QC,
         loadProbeInfo,
         differential_expression,
         barplotGO,
         pathways)
def full():
    pass

###################################################################
###################################################################
###################################################################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    to_cluster = False
    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )    

    # target = sys.argv[-1]
    # if sys.argv[-2] in ["make", "show", "plot"]:
    #     do = sys.argv[-2]
    # if "-v" in "".join(sys.argv):
    #     verbose = int([x[-1] for x in sys.argv if x[:-1] == "-v"][0])
    # else:
    #     verbose = 1
    # if do == "make":
    #     pipeline_run(target, verbose = verbose, multiprocess=8)
    # elif do == "plot":
    #     pipeline_printout_graph("pipeline.dot", "dot", target)
    # elif do == "show":
    #     pipeline_printout(sys.stdout, target, verbose = verbose)
        

