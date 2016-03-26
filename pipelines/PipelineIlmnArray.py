########################################################
########################################################
########################################################
# Functions for use with PipelineIlmnArray.py
########################################################
########################################################
########################################################

import sys, re, os
from rpy2.robjects import r as R
import collections
import numpy as np
import itertools
import sqlite3
import CGATPipelines.Pipeline as P
import random
import CGAT.Experiment as E

def normaliseAndQc( infiles, 
                    outfiles, 
                    rdir, 
                    platform = "affy", 
                    filter_above = 7, 
                    filter_over = 3,  
                    preprocess_method="gcrma",
                    preprocess = "limma", 
                    method = "quantile" ):
    '''
    python wrapper for R functions to perform
    normalisation and qc plots.
    The filter_over argument specifies the number of arrays
    that are required for a probe to be called as expressed.
    '''

    infile, ctrlfile = infiles[0], infiles[1]
    outf_postprocessed_density = outfiles[1]
    outf_postprocessed_boxplot = outfiles[2]
    outf_mean_variance = outfiles[3]

    R('''library("limma")
         library("ggplot2")
         library("affy")
         library("reshape")
         library("lumi")
         library("gcrma")
         source("%s/PipelineIlmnArray.R")''' % rdir)

    if platform == "affy":

        R('''dat <- ReadAffy()''')
        if preprocess_method == "gcrma":
            R('''dat.expression <- gcrma(dat)''')
        elif preprocess_method == "rma":
            R('''dat.expression <- rma(dat)''')
        else:
            raise ValueError("cannot perform %s preprocess" % preprocess_method)
        # build expression set
        R('''eset <- exprs(dat.expression)''')
        R('''colnames(eset) <- unlist(strsplit(colnames(eset), ".CEL"))''')
        
        # remove affy control probes
        R('''eset <- eset[grep("AFFX", rownames(eset), invert = TRUE),]''')
        # get expressed probes
        R('''expressed <- rowSums(eset > %i) >= %i''' % (filter_above, filter_over))
        R('''dat.norm.filtered <- eset[expressed, ]''')

    else:

        # read in data - maintain detection p-values for bg correction
        R('''dat.expression <- read.ilmn("%s", other.columns = "Detection")''' % infile)
        # get those that are expressed above background
        R('''expressed <- rowSums(dat.expression$other$Detection > 0.95) >= %i''' % filter_over)
        
        if preprocess == "limma":
            R('''dat.norm <- preprocessLimma("%s", "%s", method = "%s")''' % (infile, ctrlfile, method))
        elif preprocess == "lumi":
            R('''dat.norm <- preprocessLumi("%s", "%s", method = "%s")''' % (infile, ctrlfile, method))

        # remove probes expressed below background
        R('''dat.norm.filtered <- dat.norm[expressed, ]''')

    # mean vs. variance plot
    R('''stats <- getStats(dat.norm.filtered)
        pdf(file = "%s")
        mvPlot(stats)
        dev.off()''' % outf_mean_variance)

    # plots after normalisation
    R('''pdf("%s", height = 500, width = 500)
        plotDensity(dat.norm.filtered, cex.axis = 3, lwd = 3)
        dev.off()''' % outf_postprocessed_density)
      
    R('''dat.norm.filtered <- as.data.frame(dat.norm.filtered)''')
    # boxplot
    R('''gboxplot(dat.norm.filtered)
         ggsave("%s", height = ncol(dat.norm.filtered))''' % outf_postprocessed_boxplot)

    # write out expression matrix
    R('''dat.norm.filtered$ids <- rownames(dat.norm.filtered)
         write.table(dat.norm.filtered, file = "%s", sep = "\t", row.names = F)''' % outfiles[0])

########################################################
########################################################
########################################################
def buildProbeInfo(infiles, outfile):
    '''
    build probe info
    '''
    pre_infile, post_infile = infiles[0], [x for  x in infiles if x.endswith(".matrix")][0]

    outf = open(outfile, "w")
    outf.write("nprobes_prefilter\tnprobes_postfilter\n")
    ns = []

    for inf in [pre_infile, post_infile]:
        ns.append( len(open(inf).readlines()) )
    outf.write("\t".join(map(str,ns)) + "\n")
    outf.close()


########################################################
########################################################
########################################################
def buildProbe2GeneMap(infile, 
                       outfile,
                       PARAMS,
                       platform = "affy"):
    '''
    build file mapping probe id to gene id
    '''
    if platform == "affy":
        array = PARAMS.get("affy_array")
        dataset = PARAMS.get("affy_dataset")
        
        R('''library("biomaRt")''')
        R('''library("affy")''')
        R('''dat <- ReadAffy()''')

        E.info("getting probes")
        R('''probes <- featureNames(dat)''')

        E.info("getting mart")
        R('''mart <- useMart("ensembl", dataset = "%s")''' % dataset)

        # matches to hgnc symbol - this might not be appropriate for mouse data...
        E.info("mapping probes to gene")
        R('''probe2gene <- getBM(attributes = c("%s", "external_gene_name"), filters = "%s", values = probes, mart = mart)''' % (array, array))
        R('''colnames(probe2gene) <- c("probe", "gene")''')
        R('''probe2gene$gene <- toupper(probe2gene$gene)''')

        # remove probes that have no gene assignment (i.e returned "" from biomaRt) and those with 
        # multiple gene assignments - cross-hyb
        temp = P.getTempFile(".")
        E.info("writing temp file")
        R('''write.table(probe2gene, file = "%s", sep = "\t", row.names = F)''' % temp.name)
        temp.close()
        E.info("filtering probes")
        inf = open(temp.name)
        header = inf.readline()
        outf = open(outfile, "w")
        outf.write(header)
        counts = collections.defaultdict(int)
        probe2gene = {}
        for line in inf.readlines():
            data = line[:-1].split("\t")
            probe, gene = data[0], data[1]
            if gene.strip('"') == '': continue
            probe2gene[probe] = gene
            counts[probe] += 1
        for probe, count in probe2gene.iteritems():
            if count > 1:
                outf.write("%s\t%s\n" % (probe, probe2gene[probe]))
        outf.close()
        os.unlink(temp.name)
    else:
        R('''
          library(limma)
          # read in data - maintain detection p-values for bg correction
          dat <- read.ilmn(files = "%s", other.columns = "Detection")
          probe2gene <- data.frame("probe" = rownames(dat), "gene" = dat$genes$TargetID)
          write.table(probe2gene, file = "%s", row.names = F, sep = "\t")
          ''' % (infile, outfile))

########################################################
########################################################
########################################################
def clusterSamplesOnExpressionValues(infile, outfile, 
                                     dist_method = "correlation", 
                                     linkage_method = "average", 
                                     nboot = 10):
    '''
    use pvclust to cluster samples
    '''
    R('''
      library(pvclust)
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      rownames(dat) <- dat$ids
      dat <- dat[,1:ncol(dat)-1]
      clust <- pvclust(as.matrix(dat), method.hclust = "%s", method.dist = "%s", nboot = %i)
      pdf("%s")
      plot(clust)
      #, lwd = 4, cex = 4, cex.axis = 4, cex.pv = 3)
      dev.off()
      ''' % (infile, linkage_method, dist_method, int(nboot), outfile))

########################################################
########################################################
########################################################
def pcaSamplesOnExpressionValues(infile, outfile):
    '''
    use prcomp to cluster samples
    '''
    R('''
      library("ggplot2")
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      rownames(dat) <- dat$ids
      dat <- dat[,1:ncol(dat)-1]

      pc.dat <- prcomp(as.matrix(t(dat)))
      pc.dat.scores <- data.frame(pc.dat$x)

      # condition colours
      for (i in 1:nrow(pc.dat.scores)){
          name = unlist(strsplit(rownames(pc.dat.scores)[i], ".R"))
          name = name[1]
          pc.dat.scores$cond[i] <- name
      }    

      ggplot(data.frame(pc.dat.scores), aes(x = PC1, y = PC2, colour = cond)) + xlim(min(pc.dat.scores$PC1)-50, max(pc.dat.scores$PC1)+50) + geom_point(size = 6)
      ggsave(file = "%s", limitsize = F)
      ''' % (infile, outfile))

########################################################
########################################################
########################################################
def correlateSamplesOnExpressionValues(infile, outfile):
    '''
    correlate samples and draw a heatmap
    '''
    R('''
      library("pheatmap")
      library("RColorBrewer")
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      rownames(dat) <- dat$ids
      dat <- dat[,1:ncol(dat)-1]
      
      # correlations
      cors <- cor(as.matrix(dat))
      cols = brewer.pal(9, "Reds")      

      pdf("%s")
      pheatmap(cors, col = cols)
      dev.off()
      ''' % (infile, outfile))

########################################################
########################################################
########################################################
def buildDesignFile(infile, outfile, paired = False):
    '''
    build design based on conditions specified in the 
    header
    '''
    two_groups = False
    multi_group = False

    header = open(infile).readline().split("\t")
    header = [h for h in header if h != '"ids"\n' ]

    # condition associated with each sample is a function
    # of its name in the header
    conds = [re.match("(.*)(\.R[0-9]*)", x.strip('"')).groups()[0] for x in header]
    pairs = [re.match("(.*)(\.R[0-9]*)", x.strip('"')).groups()[1].split(".R")[1] for x in header]
    print pairs

    # number of distinct conditions
    nconds = len(set(conds))

    # append a list of indices that correspond
    # to each coondition
    indices = collections.OrderedDict()
    for cond in conds:
        indices[cond] = []

    i = 0
    for cond in conds:
        indices[cond].append(i)
        i += 1

    # set up design matrix and populate 
    # according to indices of conditions
    if paired:
        design = np.zeros((len(conds), nconds + 1))
        for i in range(len(conds)):
            for j in range(len(indices.keys())):
                if i in indices[indices.keys()[j]]:
                    design[i,j] = 1
                    design[i, len(design[0])-1] = pairs[i]
                else:
                    continue
                    
    else:
        design = np.zeros((len(conds), nconds))
        for i in range(len(conds)):
            for j in range(len(indices.keys())):
                if i in indices[indices.keys()[j]]:
                    design[i,j] = 1
                else:
                    continue

    # write out design file
    outf = open(outfile, "w")
    if paired:
        outf.write("\t".join(indices.keys() + ["pair"]) + "\n")
    else:
        outf.write("\t".join(indices.keys()) + "\n")
    for data in design:
        outf.write("\t".join(map(str, map(int, data))) + "\n")
    outf.close()

########################################################
########################################################
########################################################
def buildContrasts(design):
    '''
    return a list of contrasts from a design file
    '''
    conds = [h.strip() for h in open(design).readline().split("\t")]
    conds = [c for c in conds if c != "pair"]
    contrasts = []
    for cond1, cond2 in itertools.combinations(conds, 2):
        contrasts.append(cond1+"-"+cond2)
    return contrasts

########################################################
########################################################
########################################################
def runLimma(infiles, outfiles, contrasts, adjust_method = "BH", rdir = "src"):
    '''
    run limma with all conditions in the model
    '''
    data, design, probe2gene = infiles[0], infiles[1], infiles[2]
    outdir = "differential_expression.dir"
    
    # for simplicity, at the moment we gather results for each pairwise
    # comparison
    R('''
      # import libraries
      library("limma")
      source("%s/PipelineIlmnArray.R")

      # fit linear model
      dat.matrix <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      rownames(dat.matrix) <- dat.matrix$ids
      dat.matrix <- dat.matrix[,1:ncol(dat.matrix)-1]
      design <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      fit <- lmFit(dat.matrix, design)

      # make contrasts matrix
      contrasts <- pystring2rvector("%s")
      contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)

      # robust bayes fit
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      coeffs = seq(1, ncol(contrast.matrix), 1)

      # get probe2gene mapping
      probe2gene <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")

      # get results from each pairwise comparison
      for (i in coeffs){
           outname = paste("%s", paste(colnames(contrast.matrix)[i], "result", sep = "."), sep = "/")
           result <- topTable(fit2, coef = i, adjust = "%s", number = nrow(dat.matrix))
           result$ID <- rownames(result)
           result <- addGenes(result, probe2gene)
           write.table(result, file = outname, sep = "\t", row.names = F)
       }
       ''' % (rdir, data, design, contrasts, probe2gene, outdir, adjust_method))

########################################################
########################################################
########################################################
def MAPlotResults(infile, outfile, rdir = "src", p_threshold = 0.05, fc_threshold = 1):
    '''
    Produce MA plot of results
    '''
    title = os.path.basename(infile).replace(".result", "").replace("-", " vs. ")
    R('''
      source("%s/PipelineIlmnArray.R")
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      MAPlot(dat, main = "%s", fc = %f, p = %f)
      ggsave(file = "%s")
      ''' % (rdir, infile, title, fc_threshold, p_threshold, outfile))

########################################################
########################################################
########################################################
def volcanoPlotResults(infile, outfile, fc_threshold=1, p_threshold=0.05, rdir = "src"):
    '''
    Produce volcano plot of results
    '''
    title = os.path.basename(infile).replace(".result", "").replace("-", " vs. ")
    R('''
      source("%s/PipelineIlmnArray.R")
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      VolcanoPlot(dat, main = "%s", fc = %f, p = %f)
      ggsave(file = "%s")
      ''' % (rdir, infile, title, fc_threshold, p_threshold,outfile))

########################################################
########################################################
########################################################
def buildDifferentiallyExpressedGeneList(infiles, outfile, 
                                         dbh, 
                                         result_tablenames,
                                         probe2gene_table,
                                         l2fold = 1,
                                         padj = 0.05):
    '''
    annotate genes on the array as to whether they are differentially
    expressed or not. For all pairwise comparisons
    '''
    outf = open(outfile, "w")
    outf.write("contrast\tprobe_id\tgene_id\tstatus\n")
    l2fold_threshold = float(l2fold)
    p_threshold = padj
    cc = dbh.cursor()
    for table in result_tablenames:
        contrast = table
        statement = """SELECT
                       a.ID, b.gene, a.logFC, a.adj_P_Val
                       FROM %s as a, %s as b 
                       WHERE a.ID == b.probe
                       """ 
        for data in cc.execute(statement % (table, probe2gene_table)):
            probe, gene, log2fold, qvalue = data[0], data[1], data[2], data[3]
            if log2fold > l2fold and qvalue < padj:
                outf.write("\t".join(map(str, [contrast, probe, gene, 1])) + "\n")
            elif log2fold < -1*l2fold and qvalue < padj:
                outf.write("\t".join(map(str, [contrast, probe, gene, 2])) + "\n")
            else:
                outf.write("\t".join(map(str, [contrast, probe, gene, 0])) + "\n")
    outf.close()
        
########################################################
########################################################
########################################################
def heatmapDifferentiallyExpressedGeneList(dbh,
                                           diff_table,
                                           data_matrix,
                                           outfile,
                                           rdir = "src"):
    '''
    heatmap expression levels for genes that are
    differentially expressed - union of all conditions
    '''
    diff_probes = set()
    cc = dbh.cursor()
    statement = """SELECT
                   probe_id 
                   FROM %s
                   WHERE 
                   status == 1 OR status == 2"""
    for probe in cc.execute(statement % diff_table):
        diff_probes.add(probe[0])

    # for R string sub
    diff_probes = ",".join(map(str, list(diff_probes)))

    R('''
      library("pheatmap")
      source("%s/PipelineIlmnArray.R")
      library("RColorBrewer")      

      # read in and format expression matrix
      dat.matrix <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      dat.matrix <- formatDataMatrix(dat.matrix)
      
      # get the vector of probe ids
      probes = pystring2rvector("%s")

      # subset data matrix
      diff.matrix <- dat.matrix[probes,]

      # heatmap
      cols = colorRampPalette(c("blue", "white", "red"))(20)
      pdf("%s")
      pheatmap(diff.matrix, col = cols, scale = "row", fontsize = 40)
      dev.off()
      ''' % (rdir, data_matrix, diff_probes, outfile))

########################################################
########################################################
########################################################
def buildBackgroundFile(infile, outfile):
    '''
    build background sets for input into runGO.py
    '''
    outf = open(outfile, "w")
    outf.write("gene_id\n")
    for gene in set([x.split("\t")[2] for x in open(infile).readlines()]):
        outf.write(gene + "\n")
    outf.close()

########################################################
########################################################
########################################################
def buildForegroundFile(infile, outfile):
    '''
    build foreground sets for input into runGO.py
    '''
    foreground = collections.defaultdict(list)
    genes = collections.defaultdict(list)


    outf = open(outfile, "w")
    # upregulated genes
    inf = open(infile)
    header = inf.readline()
    for line in inf.readlines():
        data = line[:-1].split("\t")
        contrast, gene_id, status = data[0], data[2], data[3]
        # change status of downregulated genes
        # to zero
        if int(status) == 2:
            status = 0
        genes[contrast+"_up"].append(gene_id)
        foreground[contrast+"_up"].append(status)

    # downregulated genes
    inf = open(infile)
    header = inf.readline()
    for line in inf.readlines():
        data = line[:-1].split("\t")
        contrast, gene_id, status = data[0], data[2], data[3]
        # change status of upregulated genes
        # to zero
        if int(status) == 1:
            status = 0
        elif int(status) == 2:
            status = 1
        genes[contrast+"_down"].append(gene_id)
        foreground[contrast+"_down"].append(status)

    result = [genes.values()[0]] + foreground.values()
    outf.write("\t".join(["gene_id"] + foreground.keys()) + "\n")
    genes = set()
    for v in zip(*result):
        if v[0] in genes: continue
        genes.add(v[0])
        outf.write("\t".join(map(str, list(v))) + "\n")
    outf.close()

########################################################
########################################################
########################################################
def barplotGO(infile, outfile):
    '''
    barplot the significant pathways from GO analysis
    '''

    title = os.path.basename(infile).replace("_result.kegg.overall", "")
    R('''
      library("ggplot2")
      dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
      # sort by p-value
      dat <- dat[order(dat$fdr), ]

      # color by whether significant or not
      dat$sig <- ifelse(dat$fdr < 0.05, "yes", "no")
      dat$col <- ifelse(dat$fdr < 0.05, "blue", "red")     

      # plot the top 10
      dat <- dat[1:10,]
      plot1 <- ggplot(dat, aes(x = description, y = ratio, fill = sig, stat = "identity"))
      plot2 <- plot1 + geom_bar(stat = "identity") 
      plot3 <- plot2 + coord_flip() + theme(text=element_text(size=35, color = "black"))
      plot3 + scale_y_continuous(limits = c(0, max(dat$ratio + 1))) + scale_fill_manual(values = dat$col, labels = dat$sig)
      ggsave(file = "%s", width = 11) 
      ''' % (infile, outfile))

########################################################
########################################################
########################################################
def getPathwayGenes(goresult, gene2pathways):
    '''
    return a dictionary mapping pathways to genes
    for the top 10 pathways in pathways analysis
    '''
    go_result, gene2pathways = open(goresult), open(gene2pathways)
    h1 = go_result.readline()
    h2 = gene2pathways.readline()
    pathway2genes = collections.defaultdict(set)
    
    desc2order = []
    for line in go_result.readlines():
        data = line[:-1].split("\t")
        description, fdr = data[13], data[14]
        desc2order.append((description, fdr))
    desc2order.sort(key=lambda x: x[1])
    
    # take top 10 as in barplots
#    desc2order = desc2order[0:10]

    # map pathways to genes
    for line in gene2pathways.readlines():
        data = line[:-1].split("\t")
        description, gene = data[3], data[1]
        if description in [x[0] for x in desc2order]:
            pathway2genes[description].add(gene)

    return pathway2genes
            
########################################################
########################################################
########################################################
def getDiffGenes(dbh, contrast):
    '''
    get set of differentially expressed genes 
    '''
    diff_genes = set()
    cc = dbh.cursor()
    for data in cc.execute("""SELECT gene_id FROM differentially_expressed
                              WHERE contrast == '%s' AND status != 0""" % contrast).fetchall():
        diff_genes.add(data[0])
    return diff_genes
    
########################################################
########################################################
########################################################
def vennOverlap(dbh, outfiles):
    '''
    use limma to venn diagram the overlap between 
    two conditions
    '''
    cc = dbh.cursor()
    contrasts = set()
    for data in cc.execute("""SELECT contrast FROM differentially_expressed""").fetchall():
        contrasts.add(data[0])

    for cont1, cont2 in itertools.combinations(contrasts, 2):
        result = collections.defaultdict(list)
        temp = P.getTempFile(".")

        for data in cc.execute("""SELECT contrast, probe_id, status FROM differentially_expressed""").fetchall():
            contrast, probe, status = data[0], data[1], data[2]
            if probe in result: continue
            result["probe"].append(probe)
            if status == 2:
                status = 1
            if contrast == cont1:
                result[cont1].append(status)
            elif contrast == cont2:
                result[cont2].append(status)
        temp.write("\t".join(["probe", P.snip(cont1, "_result"), P.snip(cont2, "_result")]) + "\n")

        for data in zip(*[result["probe"], result[cont1], result[cont2]]):
            temp.write("\t".join(map(str, list(data))) + "\n")
        temp.close()
        inf = temp.name
        outf = os.path.join("differential_expression.dir", "%s_vs_%s.venn.pdf" % (cont1, cont2))
        
        R('''
        library("limma")
        dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
        venn.data <- vennCounts(dat)

        pdf("%s")
        vennDiagram(venn.data, circle.col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)), lwd = 2)
        dev.off()
        ''' % (inf, outf))

        os.unlink(inf)




########################################################
########################################################
########################################################
def buildPathwayGenes(goresult, gene2pathways, contrast, dbh, result_table, outfile):
    '''
    build matrix of fold changes for genes in pathways
    '''
    cc = dbh.cursor()
    pathway2genes = getPathwayGenes(goresult, gene2pathways)
    diff_genes = getDiffGenes(dbh, contrast)

    # just keep genes associated with pathways if they
    # are in the differentially expressed list
    for pathway, genes in pathway2genes.iteritems():
        pathway2genes[pathway] = genes.intersection(diff_genes)

    # get fold change data
    outf = open(outfile, "w")
    outf.write("pathway\tgene\tl2fold\n")

    # keep genes consistent with up/down regulation
    if goresult.find("down") != -1:
        statement = """SELECT gene, logFC FROM %s WHERE logFC < 0"""
    elif goresult.find("up") != -1:
        statement = """SELECT gene, logFC FROM %s WHERE logFC > 0"""

    # output the results
    for data in cc.execute(statement % result_table).fetchall():
        gene, l2fold = data[0], data[1]
        for pathway, g in pathway2genes.iteritems():
            if gene in g:
                outf.write("\t".join([pathway, gene, str(l2fold)]) + "\n")
    outf.close()

########################################################
########################################################
########################################################
def plotPathwayGenes(infile, outfile):
    '''
    plot the genes that are differentially expressed
    and fall into pathways
    '''
    # R will not be able to plot anything if none of the 
    # differentially expressed genes are associated
    # with a pathway. plot nothing if this is the case

    # colour of the pathways should associate with the 
    # track that they come from 

    # because the plots can get unwieldy with large 
    # gene sets, if there are more than 20 genes
    # associated with a pathway then take the top 20
    # This should be explained in the documentation

    col = random.sample(range(1,600,1), 1)[0]
    track = os.path.basename(infile).replace(".genes", "")

    if len(open(infile).readlines()) == 1:
        R('''pdf("%s")
             plot(c(0,1,2,3,4), c(0,1,2,3,4), cex = 0)
             text(2, y = 2, labels = "No genes were associated with pathways", cex = 1)
          ''' % outfile.replace(".plots", ".pdf"))
        P.touch(outfile)
    else:
        # NB. size of plot should be proportional to the
        # number of genes in the pathways
        R('''
          library("ggplot2")
          dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
          pathways <- unique(dat$pathway)
          for (p in pathways){
              toPlot <- aggregate(l2fold~gene, dat[dat$pathway == p,], mean)
              if (regexpr("/", p)[1] != -1){
                  # "/" in name not compatible with outfile names
                  p <- sub("/", "|", p)}
              outf <- paste(paste("pathways.dir/", paste("%s", p, sep = "."), sep = ""), "genes.pdf", sep = ".")
              cols <- col2rgb(%i)
              col <- rgb(cols[1], cols[2], cols[3], maxColorValue = 255)
              toPlot$col <- col
              if (nrow(toPlot) > 10){
                  toPlot <- toPlot[order(abs(toPlot$l2fold), decreasing = T),][1:10,]}
              plot1 <- ggplot(toPlot, aes(x = gene, y = l2fold, fill = col, stat = "identity")) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = toPlot$col)
              plot1 + ggtitle(p) + theme(text = element_text(size = 40, color = "black"), axis.text = element_text(colour = "Black"))
              ggsave(file = outf, width = 11, height = nrow(toPlot), limitsize = F)
          }
        ''' % (infile, track, col))
        P.touch(outfile)

        # DEPRECATED - plots were not easy to visualise genes. The code may be useful
        # in the future as an example of using geom_tile
        # R('''
        # library("ggplot2")
        # library("grid")
        # library("RColorBrewer")
        # cols = colorRampPalette(c("blue", "black", "yellow"))(20)
        # dat <- read.csv("%s", header = T, stringsAsFactors = F, sep = "\t")
        # pathways <- unique(dat$pathway)
        # list.lens <- c()
        # for (p in pathways){
        #     list.lens <- append(length(dat$pathway[dat$pathway == p]), list.lens)}
        # height <- max(list.lens)
        # height <- height + height + 10
        # plot1 <- ggplot(dat, aes(x = pathway, y = gene, fill = l2fold)) + geom_tile()
        # plot2 <- plot1 + scale_fill_gradient2(low = "yellow", mid = "black", high = "blue")
        # plot3 <- plot2 + ggtitle("%s") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        # plot3 + theme(axis.text = element_text(size = 40, color = "black"), plot.title = element_text(size = 50), legend.key.size = unit(3, "cm"), legend.text = element_text(size = 60))
        # ggsave(file = "%s", height = height, width = 20, limitsize = F)
        # ''' % (infile, title, outfile))




