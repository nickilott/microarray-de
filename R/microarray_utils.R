##########################################
##########################################
##########################################
# some generally useful functions
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
                            return(probe2gene)
                            }


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
                            return(probe2gene)
                            }



