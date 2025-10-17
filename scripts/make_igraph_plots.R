library(igraph)
library(statGraph)

args<-commandArgs(TRUE)                                                                             


el_fullpath <- file.path(getwd(), 'intermediates', paste0(args[1], ".edgelist"))                                        
meta_fullpath <- file.path(getwd(), 'intermediates', paste0(args[1], ".edgelist_meta"))
pdf_fullpath <- file.path(getwd(), 'plots', paste0(args[1], ".tumor_antigen_network.pdf"))
print(el_fullpath)
print(meta_fullpath)
p <- read_graph(el_fullpath)
node_meta <- read.csv(meta_fullpath, sep="\t", header=FALSE)


node.color <- setNames(c(node_meta$V4), c(node_meta$V1))
node.size <- setNames(c(log2(node_meta$V2 + 1)), c(node_meta$V1))
node.labels <- setNames(c(node_meta$V3), c(node_meta$V1))

legend_df <- unique(node_meta[c("V4", "V5")])
legend_df <- legend_df[complete.cases(legend_df),]

print(legend_df)

entropy <- graph.entropy(p)

print(entropy)

pdf(pdf_fullpath)
plot(p, xlim=c(-1, 1.75), edge.arrow.size = 0, vertex.size=node.size, vertex.color=node.color, vertex.label=node.labels, vertex.label.cex = 0.001)
legend(x=1, y=1, legend=c(legend_df$V5), fill=c(legend_df$V4), cex=0.65, bty="n")
title(paste0("pMHC Network Among Tumor Cells\n", args[1], "\nEntropy: ", round(entropy, 2)), line=-3)
dev.off()
