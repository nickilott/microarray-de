##########################################
##########################################
##########################################
# some generally useful functions
##########################################
##########################################
##########################################

# load packages
library("ggplot2")
library("grid")
library("plyr")
library("ggrepel")

##########################################
##########################################
##########################################

readData <- function(tsv, header=TRUE){
	 return (read.csv(tsv,
	                 header=header,
			 stringsAsFactors=F,
			 sep="\t"))
			 }

##########################################
##########################################
##########################################

filterNonUniqueGenes <- function(probe2gene){
                            genes = c()
                            probes = c()
                                     
                            for (i in 1:nrow(probe2gene)){
                                 gene <- as.character(probe2gene[i,1])
                                 if (length(gene) == 0){next}
                                 probe <- as.character(probe2gene[i,2])
                                 if (!(gene %in% genes)){
                                    genes <- append(genes, gene)
                                    probes <- append(probes, probe)
                                 }else
                                 {
                                 next
                                 }
                             }
                            probe2gene <- probe2gene[probes,]
			    probe2gene[,1] <- as.character(probe2gene[,1])
			    probe2gene[,2] <- as.character(probe2gene[,2])	
			    return(na.omit(probe2gene))
                            }


##########################################
##########################################
##########################################

filterNonUniqueProbes <- function(probe2gene){
                            genes = c()
                            probes = c()
                            index = 0
                            keep = c()
                            for (i in 1:nrow(probe2gene)){
                                 index = index + 1
                                 gene <- probe2gene[i,1]
                                 probe <- as.character(probe2gene[i,2])
                                 if (!(probe %in% probes)){
                                    genes <- append(genes, gene)
                                    probes <- append(probes, probe)
                                    keep <- append(keep, index) 
                                 }else
                                 {
                                 next
                                 }
                             }
                            probe2gene <- probe2gene[keep,]
			    probe2gene[,1] <- as.character(probe2gene[,1])
			    probe2gene[,2] <- as.character(probe2gene[,2])	
			    return(na.omit(probe2gene))
}

##########################################
##########################################
##########################################

runPCA <- function(matrix.file){

       # read in and format data
       dat <- readData(matrix.file)
       rownames(dat) <- dat$ids
       dat <- dat[,1:ncol(dat)-1]

       # run PCA
       pc.dat <- prcomp(as.matrix(t(dat)))

       # return pca object
       return (pc.dat)
}


##########################################
##########################################
##########################################


buildPCAScores <- function(pc.dat, outfile){
	       # return pca scores
	       scores <- data.frame(pc.dat$x)
	       scores$sample <- rownames(scores)
	       scores <- scores[,c("sample", colnames(scores)[1:ncol(scores)-1])]
	       write.table(scores,
	                   file=outfile,
			   sep="\t",
			   row.names=F,
			   quote=F)
}

##########################################
##########################################
##########################################

buildPCAVarianceExplained <- function(pc.dat, outfile){
			  # get variance explained from
			  # pca
			  ve <- data.frame(summary(pc.dat)$importance)
			  ve <- data.frame(t(ve[2,]))
			  ve$PC <- rownames(ve)
			  ve <- ve[,c("PC", "Proportion.of.Variance")]
			  write.table(ve,
				      file=outfile,
				      sep="\t",
				      row.names=F,
				      quote=F)
}				      


##########################################
##########################################
##########################################

plotPCA <- function(pc.dat,
                    pcs=c("PC1", "PC2"),
		    outfile="outfile_pca.pdf"){

      # get scores
      pc.dat.scores <- data.frame(pc.dat$x)

      # condition colours
      for (i in 1:nrow(pc.dat.scores)){
	  name = unlist(strsplit(rownames(pc.dat.scores)[i], ".R"))
	  name = name[1]
          pc.dat.scores$cond[i] <- name
      }

      # get variance explained
      imps <- data.frame(summary(pc.dat)$importance)
      pc1 <- pcs[1]
      pc2 <- pcs[2]
      imps <- c(imps[,pc1][2], imps[,pc2][2])

      minimum = min(append(pc.dat.scores[,pc1], pc.dat.scores[,pc2]))
      maximum = max(append(pc.dat.scores[,pc1], pc.dat.scores[,pc2]))

      plot1 <- ggplot(data.frame(pc.dat.scores), aes(x=get(pcs[1]), y=get(pcs[2]), colour = cond, fill=cond))
      plot2 <- plot1  + geom_point(shape=18, size = 6) 
      plot2 + xlab(as.character(imps[1])) + ylab(as.character(imps[2])) + theme_bw()
      ggsave(file = outfile, limitsize = F)
#+ stat_ellipse(type="t", linetype="dashed", geom="polygon", level=0.95, alpha=0.2)
}

##########################################
##########################################
##########################################

plotLoadings <- function(pc.dat,
	                 probe2gene,
	     		 pcs=c("PC1", "PC2"),
	     		 outfile="outfile_loadings.pdf",
	     		 top=5,
	     		 x.lim=c(-0.15,0.15),
	     		 y.lim=c(-0.15,0.15)){

               # read probe2gene map
	       probe2gene <- readData(probe2gene)
	       rownames(probe2gene) <- probe2gene$probe

	       # get loadings
  	       loads <- data.frame(pc.dat$rotation)
 	       loads$id <- rownames(loads)
	       	       
	       loads$gene <- probe2gene[loads$id,]$gene

	       # remove non-unique genes
	       to.filter <- data.frame(cbind(loads$gene,loads$id))
	       to.filter[,1] <- as.character(to.filter[,1])
	       to.filter[,2] <- as.character(to.filter[,2])	
	       rownames(to.filter) <- to.filter[,2]

	       filtered <- filterNonUniqueGenes(to.filter)
	       print (dim(filtered))

	       loads <- loads[filtered[,2],]

	       # set upi where arrows are drawn from 
               loads$x <- 0
	       loads$y <- 0

	       # get the top specified for labeling
	       toppc1 <- rownames(loads[order(-loads[,pcs[1]]),])[1:top]
	       bottompc1 <- rownames(loads[order(loads[,pcs[1]]),])[1:top]

	       toppc2 <- rownames(loads[order(-loads[,pcs[2]]),])[1:top]
	       bottompc2 <- rownames(loads[order(loads[,pcs[2]]),])[1:top]

	       to.keep <- c(toppc1, bottompc1, toppc2, bottompc2)
	       
	       # subset the data
	       loads <- loads[to.keep,]

 	       plot1 <- ggplot(loads, aes(x=x, y=y, xend=get(pcs[1]), yend=get(pcs[2])))
 	       plot2 <- plot1 + geom_segment(arrow=arrow(length=unit(0.5, "cm"), type="closed"))
	       plot3 <- plot2 + xlim(x.lim[1], x.lim[2]) + ylim(y.lim[1], y.lim[2])
	       plot3 + geom_text(data=loads, aes(x=get(pcs[1]), y=get(pcs[2]), label=gene), size=8) + guides(colour=FALSE)
 	       ggsave(outfile)
	       
}

##########################################
##########################################
##########################################

summariseData <- function(dat, cond1, cond2){
	         df <- ddply(dat, c("probe", "cond"), summarise, average=mean(value))
	         df <- data.frame(cbind(df[df$cond == cond1,],
	                                   df[df$cond == cond2,]))
		 return(df)
}

##########################################
##########################################
##########################################

annotateGenes <- function(df, anno){
    		# annotate genes
        	df$status <- ifelse(df$probe %in% anno$probe, "annotated", "not changed")
		for (p in anno$probe){
	            df$annotation[df$probe == p] <- anno$gene[anno$probe == p]
		    }
		return (df)

}

##########################################
##########################################
##########################################

getClusterColours <- function(clusters){
	             # get colours associated with cluster
	             clusters$colour <- rainbow(max(clusters$cluster),
		                                s=1,
						v=0.75)[clusters$cluster]
		     return(clusters)
}

##########################################
##########################################
##########################################

addColours <- function(df, clustersWithColours){
	      # add colours
	      df$probe <- as.character(df$probe)
	      df$colour <- clustersWithColours[df$probe,]$colour
	      df$colour <- ifelse(is.na(df$colour), "grey", df$colour)
	      return(df)
}

##########################################
##########################################
##########################################

addClusters <- function(df, clusters){
		df$cluster <- clusters[df$probe,]$cluster
		df$cluster <- ifelse(is.na(df$cluster),
		                     "not changed",
				      df$cluster)
		return(df)
}

##########################################
##########################################
##########################################

scatterplotAbundances <- function(matrix.file,
                                  annotation,
			          cond1,
				  cond2,
				  outfile){

 		         library(reshape)
 		         library(plyr)

		         # read annotations
			 anno <- readData(annotation)

 			 # read matrix
        		 dat <- readData(matrix.file)
			 rownames(dat) <- dat$ids
			 dat <- dat[,1:ncol(dat)-1]

			 # get conditions
			 conds <- unlist(strsplit(colnames(dat), ".R[0-9]"))
			 dat <- dat[,which(conds == cond1 | conds == cond2)]

			 # reformat
			 dat$probe <- rownames(dat)
			 dat <- melt(dat)

			 dat$cond <- unlist(strsplit(as.character(dat$variable), ".R[0-9]"))

			 # summarise data per probe/condition
			 dat.summarised <- summariseData(dat, cond1, cond2)

       			 # annotate genes
			 dat.summarised <- annotateGenes(dat.summarised, anno)

			 # scatterplot the results
                         plot1 <- ggplot(dat.summarised, aes(x=average, y=average.1, colour=status, label=annotation))
	                 plot2 <- plot1 + geom_point(shape=18) + xlab(cond1) + ylab(cond2)
			 plot3 <- plot2 + geom_abline(slope=1, intercept=1, linetype="dashed")
			 plot4 <- plot3 + geom_abline(slope=1, intercept=-1, linetype="dashed") 
			 plot5 <- plot4 + geom_abline(slope=1, intercept=0)
			 plot5 + geom_text_repel(colour="black") + scale_colour_manual(values=c("red", "grey"))
		         ggsave(outfile)
}

##########################################
##########################################
##########################################

scatterplotAbundancesClusters <- function(matrix.file,
                                  annotation,
				  clusters,
			          cond1,
				  cond2,
				  outfile){
				  
 		         library(reshape)
 		         library(plyr)

		         # read annotations
			 anno <- readData(annotation)

 			 # read matrix
        		 dat <- readData(matrix.file)
			 rownames(dat) <- dat$ids
			 dat <- dat[,1:ncol(dat)-1]

                         # read clusters
			 clusters <- readData(clusters)
			 rownames(clusters) <- clusters$probe

			 # get colours associated with cluster
			 clusters <- getClusterColours(clusters)

			 # get conditions
			 conds <- unlist(strsplit(colnames(dat), ".R[0-9]"))
			 dat <- dat[,which(conds == cond1 | conds == cond2)]

			 # reformat
			 dat$probe <- rownames(dat)
			 dat <- melt(dat)

			 dat$cond <- unlist(strsplit(as.character(dat$variable), ".R[0-9]"))

			 # summarise data per probe/condition
			 dat.summarised <- summariseData(dat, cond1, cond2)

    			 # annotate genes
        		 dat.summarised <- annotateGenes(dat.summarised, anno)

			 # add colours
			 dat.summarised <- addColours(dat.summarised, clusters)

			 # add clusters
			 dat.summarised <- addClusters(dat.summarised, clusters)
 
			 # order for colours to work
		 	 dat.summarised <- dat.summarised[order(dat.summarised$cluster),]

			 # get size - big for in cluster
			 dat.summarised$size <- ifelse(dat.summarised$cluster == "not changed", 1, 4)

			 # scatterplot the results
			 col <- factor(dat.summarised$cluster, levels = sort(unique(dat.summarised$cluster)))
                         plot1 <- ggplot(dat.summarised, aes(x=average, y=average.1, colour=col, label=annotation, stat="identity", size=size))
	                 plot2 <- plot1 + geom_point(shape=18) + xlab(cond1) + ylab(cond2)
			 plot3 <- plot2 + geom_abline(slope=1, intercept=1, linetype="dashed")
			 plot4 <- plot3 + geom_abline(slope=1, intercept=-1, linetype="dashed") 
			 plot5 <- plot4 + geom_abline(slope=1, intercept=0)
			 plot5 + geom_text_repel(colour="black") + scale_colour_manual(values=unique(dat.summarised$colour))
		         ggsave(outfile)
}

##########################################
##########################################
##########################################

scatterplotFoldChanges <- function(results1,
		                  results2,	
		        	  annotation,
		       		  clusters,
		       		  outfile="scatterplot_fold.pdf"){

		       # read and format
                       results1 <- readData(results1)
		       results2 <- readData(results2)
		       rownames(results1) <- results1$ID
		       rownames(results2) <- results2$ID

		       # combine into single dataframe
                       results1 <- results1[rownames(results2),]

		       results <- data.frame(cbind(results1, results2))
		       results$probe <- results$ID.1

		       # read annotations
		       anno <- readData(annotation)

                       # read clusters
		       clusters <- readData(clusters)
		       rownames(clusters) <- clusters$probe

    		       # annotate genes
		       results <- annotateGenes(results, anno)
		       rownames(results) <- results$probe

		       # add colours
		       clusters <- getClusterColours(clusters)
		       results <- addColours(results, clusters)	

		       # assign clusters
		       results <- addClusters(results, clusters)

		       # order for colours to work
		       results <- results[order(results$cluster),]

		       # get size - big for in cluster
		       results$size <- ifelse(results$cluster == "not changed", 1, 2)

		       # plot
		       maximum <- max(append(results$logFC, results$logFC.1))
		       minimum <- min(append(results$logFC, results$logFC.1))

		       col <- factor(results$cluster, levels = sort(unique(results$cluster)))	
		       plot1 <- ggplot(results, aes(x=logFC, y=logFC.1, colour=col, label=annotation, stat="identity", size=size))
		       plot2 <- plot1 + geom_point(shape=18) + scale_size_area(max_size=3)
		       plot3 <- plot2 + geom_text_repel(colour="black", size=5) + scale_colour_manual(values=unique(results$colour))
		       plot4 <- plot3 + xlim(minimum, maximum) + ylim(minimum, maximum)
		       plot5 <- plot4 + geom_hline(yintercept=c(-1,1), linetype="dashed") + geom_vline(xintercept=c(-1,1), linetype="dashed")
		       plot5 + theme_bw() + ylim(c(min(results$logFC.1), max(results$logFC.1)))
		       ggsave(outfile)
}

##########################################
##########################################
##########################################


plotPathways <- function(pathway_results, colour="red", outfile="pathways.pdf"){
	     dat <- readData(pathway_results)
	     dat <- dat[dat$code == "+",]
	     dat <- dat[order(dat$ratio),]
	     plot1 <- ggplot(dat, aes(x=factor(description, levels=dat$description), y=ratio, stat="identity"))
	     plot2 <- plot1 + geom_bar(stat="identity", fill=colour)
	     plot2 + scale_fill_manual(values=c(colour)) + coord_flip() + theme_bw()
	     ggsave(outfile)
}


summariseByDesign <- function(matrix, design, probe){

 		  # summarise probe expression
		  # levels by conditions given
		  # in the design file

		  # read and format matrix
		  mat <- readData(matrix)
		  rownames(mat) <- mat$ids
		  mat <- mat[,grep("ids", colnames(mat), invert=T)]

		  # read design file
		  design <- readData(design)

		  # subset mat
		  mat <- mat[,design[,1]]

		  # summarise data by design
                  res <- data.frame(cbind(design[,2], design[,3], t(mat[probe,])))
		  colnames(res) <- c("strain", "stimulation", "exprs")
        	  res$exprs <- as.numeric(as.character(res$exprs))
		  res$strain <- as.character(res$strain)
		  res$stimulation <- as.character(res$stimulation)
		  res.sum <- ddply(res, 
                         c("strain", "stimulation"),
                         summarize,
                         "mean.exprs" = mean(exprs),
                         "sd" = sd(exprs),
			 "n" = length(exprs),
			 "se" = sd/sqrt(n))
        res.sum <- res.sum[order(res.sum$stimulation, decreasing = T),]
	results <- list("summarised" = res.sum, "raw" = res)
	return (results)
}

lineplotSummarised <- function(summarised, outfile, colours=c()){
          	   # takes as input summarisedByDesign
		   # output and plots by condition

		   res.sum <- summarised$summarised
		   res <- summarised$raw
		   plot1 <- ggplot(res.sum,
		                   aes(x=factor(stimulation, levels = unique(stimulation)),
				   y=mean.exprs,
				   group=strain,
				   colour = strain,
				   stat = "identity"))
		   plot2 <- plot1 + geom_line()
		   plot3 <- plot2 + geom_errorbar(aes(ymax=mean.exprs + se,
		   	    	    		      ymin=mean.exprs - se),
						      width = 0.2)
        	   plot4 <- plot3 + scale_colour_manual(values = c("blue", "red"))
		   plot4 + geom_jitter(data=res,
		                       aes(x=factor(stimulation, levels=unique(stimulation)),
				       y=exprs,
				       group=strain),
				       height=0.2,
				       width=0.2,
				       pch=18) + theme_bw()
	ggsave(outfile)
	}

boxplotSummarised <- function(summarised, outfile){
		   res.sum <- summarised$summarised
		   res <- summarised$raw
		   plot1 <- ggplot(res,
		   	           aes(x=factor(strain, levels=unique(strain)),
				   y=exprs))
		   plot2 <- plot1 + geom_boxplot()
		   plot3 <- plot2 +  geom_jitter(data=res,
		                       aes(x=factor(strain, levels=unique(strain)),
				       y=exprs),
				       height=0.2,
				       width=0.2,
				       shape=18,
				       size=3)
	           plot4 <- plot3 + facet_grid(~stimulation)
		   plot4 + theme_bw()
		   ggsave(outfile, height=7, width=14)
		   }