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


readData <- function(tsv, header=TRUE){
	 return (read.csv(tsv,
	                 header=header,
			 stringsAsFactors=F,
			 sep="\t"))
			 }

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


plotPCA <- function(pc.dat,
                    pcs=c("PC1", "PC2"),
		    outfile="outfile_pca.pdf",
		    x.lim=c(-25,25),
		    y.lim=c(-25,25)){

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

      plot1 <- ggplot(data.frame(pc.dat.scores), aes(x=get(pcs[1]), y=get(pcs[2]), colour = cond))
      plot2 <- plot1 + xlim(x.lim[1], x.lim[2]) + ylim(y.lim[1],y.lim[2]) + geom_point(shape=18, size = 6)
      plot2 + xlab(as.character(imps[1])) + ylab(as.character(imps[2]))
      ggsave(file = outfile, limitsize = F)
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
			 dat.summarised <- ddply(dat, c("probe", "cond"), summarise, average=mean(value))
	    		 dat.summarised <- data.frame(cbind(dat.summarised[dat.summarised$cond == cond1,],
	                                                dat.summarised[dat.summarised$cond == cond2,]))

    			# annotate genes
        		dat.summarised$status <- ifelse(dat.summarised$probe %in% anno$probe, "annotated", "not changed")

			for (p in anno$probe){
	                     dat.summarised$annotation[dat.summarised$probe == p] <- anno$gene[anno$probe == p]}

			# scatterplot the results
                        plot1 <- ggplot(dat.summarised, aes(x=average, y=average.1, colour=status, label=annotation))
	                plot2 <- plot1 + geom_point(shape=18) + xlab(cond1) + ylab(cond2)
			plot3 <- plot2 + geom_abline(slope=1, intercept=1, linetype="dashed")
			plot4 <- plot3 + geom_abline(slope=1, intercept=-1, linetype="dashed") 
			plot5 <- plot4 + geom_abline(slope=1, intercept=0)
			plot5 + geom_text() + scale_colour_manual(values=c("red", "grey"))
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
			 clusters$colour <- rainbow(max(clusters$cluster))[clusters$cluster]

			 # get conditions
			 conds <- unlist(strsplit(colnames(dat), ".R[0-9]"))
			 dat <- dat[,which(conds == cond1 | conds == cond2)]

			 # reformat
			 dat$probe <- rownames(dat)
			 dat <- melt(dat)

			 dat$cond <- unlist(strsplit(as.character(dat$variable), ".R[0-9]"))

			 # summarise data per probe/condition
			 dat.summarised <- ddply(dat, c("probe", "cond"), summarise, average=mean(value))
	    		 dat.summarised <- data.frame(cbind(dat.summarised[dat.summarised$cond == cond1,],
	                                                dat.summarised[dat.summarised$cond == cond2,]))

    			# annotate genes
        		dat.summarised$status <- ifelse(dat.summarised$probe %in% anno$probe, "annotated", "not changed")

			for (p in anno$probe){
	                     dat.summarised$annotation[dat.summarised$probe == p] <- anno$gene[anno$probe == p]}
			rownames(dat.summarised) <- dat.summarised$probe


			# add colours
			dat.summarised$colour <- clusters[dat.summarised$probe,]$colour
			dat.summarised$colour <- ifelse(is.na(dat.summarised$colour), "grey", dat.summarised$colour)
			dat.summarised$cluster <- clusters[dat.summarised$probe,]$cluster
			dat.summarised$cluster <- ifelse(is.na(dat.summarised$cluster), "not changed", dat.summarised$cluster)

			# order for colours to work
			dat.summarised <- dat.summarised[order(dat.summarised$cluster),]

			# scatterplot the results
			col <- factor(dat.summarised$cluster, levels = sort(unique(dat.summarised$cluster)))
                        plot1 <- ggplot(dat.summarised, aes(x=average, y=average.1, colour=col, label=annotation, stat="identity"))
	                plot2 <- plot1 + geom_point(shape=18) + xlab(cond1) + ylab(cond2)
			plot3 <- plot2 + geom_abline(slope=1, intercept=1, linetype="dashed")
			plot4 <- plot3 + geom_abline(slope=1, intercept=-1, linetype="dashed") 
			plot5 <- plot4 + geom_abline(slope=1, intercept=0)
			plot5 + geom_text(colour="black") + scale_colour_manual(values=unique(dat.summarised$colour))
		        ggsave(outfile)
		    }