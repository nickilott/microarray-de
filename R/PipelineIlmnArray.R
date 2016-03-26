########################################################
########################################################
########################################################
# functions for use with pipeline_ilmnarray.py
########################################################
########################################################
########################################################

library("ggplot2")
library("limma")
library("reshape")
library("lumi")
library("affy")

########################################################
# preprocessing dependent on package specified
########################################################

# limma
preprocessLimma <- function(infile, ctrlfile, method = "quantile"){
     
     dat <- read.ilmn(infile, ctrlfiles=ctrlfile, other.columns = "Detection")
     dat.bgc <- nec(dat)
     # remove negative controls before normalization
     dat.bgc <- dat.bgc[dat.bgc$genes$Status == "regular",]
     dat.norm <- normalizeBetweenArrays(dat.bgc, method = method)
     return (dat.norm$E)
}

# lumi
preprocessLumi <- function(infile, ctrlfile, method = "vsn"){
    
     library("vsn") 
     dat <- lumiR(infile)
     controlData(dat) <- ctrlfile
     dat.bgc <- lumiB(dat, method = "bgAdjust")

     # VST normalization
     dat.transformed <- lumiT(dat.bgc)
     dat.norm <- lumiN(dat.transformed, method = method)
     return (exprs(dat.norm))
}

########################################################
# ggplot boxplot from ilmn object
########################################################

gboxplot <- function(dat, outfile){
           dat <- melt(dat)
           p <- ggplot(dat, aes(x = variable, y = value, fill = variable))
           p + geom_boxplot(position = "dodge") + coord_flip() 
	   
}

########################################################      
# ranked averge and variance for variance plot
########################################################

getStats <- function(dat.matrix){
            ave <- apply(dat.matrix,1,mean)
            sd <- apply(dat.matrix,1,sd)
            rank <- rank(ave)
            result <- data.frame("rank" = rank, "sd" = sd)
            return(result)
}

########################################################
# plot mean vs. variance
########################################################

mvPlot <- function(stats.matrix){
          smoothScatter(stats.matrix[, "rank"], stats.matrix[, "sd"], main = "mean vs. sd", cex.axis = 2)
          lines(lowess(stats.matrix[, "rank"], stats.matrix[, "sd"]), col = "red", lty = 2, lwd = 2)
}

########################################################
# add gene names to limma results table
########################################################

addGenes <- function(result, probe2gene){
         result <- merge(result, probe2gene, by.x = "ID", by.y = "probe", all.x = T)
         rownames(result) <- result$ID
         return(result)
}

########################################################
# MA plot function
########################################################

MAPlot <- function(result, main = "", fc = 1, p = 0.05){
          result$sig <- ifelse(result$adj.P.Val < p & abs(result$logFC) > fc , 1, 0)
          plot1 <- ggplot(result, aes(x = AveExpr, y = logFC, colour = factor(sig)))
          plot2 <- plot1 + geom_point() + scale_colour_manual(values = c("black", "red")) 
          plot3 <- plot2 + ggtitle(main) + xlab("Expression level") + ylab("log2(fold change)")
          plot4 <- plot3 + geom_hline(yintercept = c(-1, 1), linetype = "dashed")
	  plot4 
}

########################################################
# Volcano plot function
########################################################

VolcanoPlot <- function(result, main = "", fc=1, p=0.05){
          result$sig <- ifelse(result$adj.P.Val < p & abs(result$logFC) > fc, 1, 0)
          plot1 <- ggplot(result, aes(x = logFC, y = -log10(P.Value), colour = factor(sig)))
          plot2 <- plot1 + geom_point() + scale_colour_manual(values = c("black", "red")) 
          plot3 <- plot2 + ggtitle(main) + xlab("log2(Fold change)") + ylab("-log10(p-value)")
	  plot4 <- plot3 + geom_vline(xintercept = c(-1, 1), linetype = "dashed") 
	  plot5 <- plot4 + geom_hline(yintercept = 2, linetype = "dashed", colour = "red")
	  plot6 <- plot5 + geom_hline(yintercept = 5, linetype = "dashed", colour = "blue")
	  plot7 <- plot6 + geom_text(data = NULL, x = max(result$logFC) - 0.5, y = 2.4, label = "p = 0.01", colour = "red" )
	  plot8 <- plot7 + geom_text(data = NULL, x = max(result$logFC) - 0.5, y = 5.4, label = "p = 1E-05", colour = "blue" )
	  plot8
}

########################################################
# format data matrix that contains ids - make ids 
# rownames and remove ids column
########################################################

formatDataMatrix <- function(dat.matrix){
		 rownames(dat.matrix) <- dat.matrix$ids
                 dat.matrix <- dat.matrix[, 1:ncol(dat.matrix) - 1]
 		 return (dat.matrix)
}

########################################################
# convert python string of the form "a,b,c,d" into 
# R vector of the form c("a","b","c","d")
########################################################

pystring2rvector <- function(string){
              	    rvector <- c(unlist(strsplit(string, ",")))
		    return(rvector)
}

########################################################
# Barplot expression values for a set of probes
########################################################

barplotProbeset <- function(dat.matrix, probe2gene, probes){

		    library(reshape)        	    

                    # gene to probe map
                    info <- probe2gene[probes, ]
                   
                    # probe expression values
                    exprs <- dat.matrix[probes, ]

                    # conditions
		    conds <- unlist(lapply(colnames(exprs), function(x){ gsub(".R[0-9]", "", x) }))

		    # take means and standard errors for each probe
		    means <- data.frame(apply(exprs, 1, function(x) { tapply(x, conds, mean) }))

      		    means$cond <- rownames(means)
		    means <- melt(means)
 		    colnames(means) <- c("condition", "probe", "mean")
		    means$probe <- unlist(lapply(means$probe, function(x){ gsub("X", "", x) }))
		    means$gene <- probe2gene[means$probe,][,2]

		    # standard errors
		    se <- function(x){
                          stdev <- sd(x)
                          return (stdev / sqrt(length(x)))}
		    
                    ses <- data.frame(apply(exprs, 1, function(x) { tapply(x, conds, se) }))

      		    ses$cond <- rownames(ses)
		    ses <- melt(ses)
 		    colnames(ses) <- c("condition", "probe", "se")

		    # combine data
		    res <- data.frame(cbind(means, ses))
		    
		    # take out NAs
		    res <- na.omit(res)
		    res$probe <- paste(res$gene, res$probe, sep = "_")

		    # plot the data
		    limits = aes(ymin = mean - se, ymax = mean + se)
		    plot1 <- ggplot(res, aes(x = probe, y = mean, fill = gene, colour = condition))
                    plot2 <- plot1 + geom_bar(position = position_dodge(), stat = "identity")
		    plot3 <- plot2 + geom_errorbar(limits, position = position_dodge(0.9), width = 0.30)
                    plot3 + coord_flip() + scale_colour_manual(values = c("black", "red"))}








