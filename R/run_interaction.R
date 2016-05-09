##############################################################
##############################################################
##############################################################
# This script is designed for microarray data where the 
# comparison of interest is for example bewteen WT and mutant
# strains that are stimulated and unstimulated. I.e. we are
# interested in the differential response of the different
# strains to stimulation
##############################################################
##############################################################
##############################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("plyr"))

# make options list
option_list <- list(
               make_option(c("-m", "--data-matrix"),
                           help="input normalised matrix"),
               make_option(c("-p", "--probe2gene"),
                           help="probe2gene mapping file"),
               make_option(c("-d", "--design"),
                           help="input design file"),
               make_option(c("-c", "--contrasts"), type="character",
                           help="specify contrasts of interest"),
               make_option(c("-n", "--probe-column"), type="character",
                           help="specify name of probe column"),
	       make_option(c("--plot-interactions"), action="store_true",
                           help="plot interactions?")
               )

##############################
# get command line options
##############################

opt <- parse_args(OptionParser(option_list=option_list))

##############################
# function for reading data in
##############################

readData <- function(infile){
	 dat <- read.csv(infile,
                         header = T, 
                         stringsAsFactors = F,
                         sep = "\t")
	return (dat)}

##############################
# read in data
##############################

mat <- readData(opt$`data-matrix`)

# user defined probe name column
probe.colname <- opt$`probe-column`
r <- mat[,probe.colname]
mat <- mat[,grep(probe.colname, colnames(mat), invert=T)]
rownames(mat) <- r

# assume probe is called "probe"
# in file
probe2gene <- readData(opt$`probe2gene`)
rownames(probe2gene) <- probe2gene$probe

# read design file
design <- readData(opt$`design`)

################################
# if numeric this is a problem
################################
#design$track <- paste("X", design$track, sep="")
##################################
##################################

# build targets
ts <- paste(design[,2], design[,3], sep = ".")
ts <- factor(ts, levels=unique(ts))

# build design
des <- model.matrix(~0+ts)
colnames(des) <- levels(ts)

# rows in design match must match cols in matrix file
# in design
tracks <- design$track
print (intersect(colnames(mat), tracks))

mat <- mat[,unlist(tracks)]


# fit model
fit1 <- lmFit(mat, des)

# get contrasts
conts <- opt$`contrasts`
conts <- unlist(strsplit(conts, ","))

cont.matrix <- makeContrasts(conts[1],
	                     conts[2],
			     conts[3],
	                     levels = des)


# fit contrasts 
fit2 <- contrasts.fit(fit1, cont.matrix)
fit2 <- eBayes(fit2)

###########################
# get the results for each 
# contrast
###########################

for (c in seq(1,length(conts),1)){
    outname <- paste(conts[c], "result", sep = ".")
    result <- topTable(fit2, coef=c, number = nrow(mat))
    result$ID <- rownames(result)
    result$gene <- probe2gene[result$ID,]$gene
    write.table(result, file = outname, sep = "\t", row.names = F)
    }

################
# plot results
################

if (opt$`plot-interactions` == TRUE){
    result <- topTable(fit2, coef = 3, number = nrow(mat))
    result$ID <- rownames(result)
    result$gene <- probe2gene[result$ID,]$gene
    diff.probes <- result$ID[result$adj.P.Val < 0.05]
    mat2 <- mat[diff.probes,]
    for (i in 1:nrow(mat2)){
        res <- data.frame(cbind(design[,2], design[,3], t(mat2[i,])))
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
        outname = paste(diff.probes[i], probe2gene[probe2gene$probe == diff.probes[i],]$gene, sep = "_")
	outname = paste(outname, "interaction.pdf", sep = "_")
	plot1 <- ggplot(res.sum, aes(x = factor(stimulation, levels = stimulation), y = mean.exprs, group = strain, colour = strain, stat = "identity"))
	plot2 <- plot1 + geom_line()
	plot3 <- plot2 + geom_errorbar(aes(ymax=mean.exprs + se, ymin=mean.exprs - se), width = 0.2)
        plot3 + scale_colour_manual(values = c("blue", "red")) + geom_jitter(data = res, aes(x=factor(stimulation, levels=stimulation), y=exprs, group=strain))
	ggsave(outname)}
	warnings()
}
