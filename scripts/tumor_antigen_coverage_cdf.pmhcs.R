library("ggplot2")
library("dplyr")
library("tidyr")

args<-commandArgs(TRUE)

f <- function(m) {
  n <- nrow(m)
  df <- data.frame(pmhc = character(n), cumulative_coverage_cell_count = integer(n))
  cs <- 0L
  
  for (i in 1:n) {
    rs <- rowSums(m)
    j <- which.max(rs)
    cs <- cs + rs[j]
    df[i, 1] <- rownames(m)[j]
    df[i, 2] <- cs
    
    if (rs[j] == ncol(m)) {
      df[(i + 1L):n, 1] <- rownames(m)[-j]
      df[(i + 1L):n, 2] <- cs
      break
    }
    
    m <- m[,m[j,] == 0, drop = FALSE][-j,,drop = FALSE]
  }
  
  df
}

find_elbow_derivative <- function(y) {
  # Calculate first and second derivatives
  first_deriv <- diff(y)
  second_deriv <- diff(first_deriv)
  
  # Find maximum of second derivative (for decreasing curves)
  # or minimum for increasing curves
  elbow_idx <- which.max(abs(second_deriv)) + 1  # +1 because diff reduces length
  
  return(y[elbow_idx])
}

counts_log <- read.delim(args[1], sep='\t', header=TRUE, row.names=1, encoding = "UTF-8", fileEncoding = "UTF-8")
print(head(counts_log))

df_proportions <- data.frame(
  pmhc = rownames(t(counts_log)),
  proportion_cells_expressing = rowMeans(t(counts_log), na.rm = TRUE)*100
)

#print(head(t(counts_log)))
pmhc_select <- f(t(counts_log))
pmhc_select$cumulative_coverage_percent <- pmhc_select$cumulative_coverage_cell_count/nrow(counts_log) * 100
pmhc_select <- merge(pmhc_select, df_proportions, by="pmhc")
pmhc_select$antigen_source <- ifelse(startsWith(pmhc_select$pmhc, "ERV."), "ERV",
                              ifelse(startsWith(pmhc_select$pmhc, "chr"), "SNV", "CTA"))

print(max(pmhc_select$cumulative_coverage_percent))
print(find_elbow_derivative(pmhc_select$cumulative_coverage_percent))
df_sorted <- pmhc_select[order(pmhc_select$cumulative_coverage_percent), ]
#first_row <- which(df_sorted$CPrec >= (find_elbow_derivative(pmhc_select$CPrec)))[1]
#first_row <- which(df_sorted$CPrec >= (max(pmhc_select$CPrec) - (max(pmhc_select$CPrec) * .1)))[1]
#pmhc_select <- df_sorted[1:first_row, ]

pmhc_select_out <- pmhc_select %>% arrange(cumulative_coverage_percent)
pmhc_select_out$pmhc <- gsub("newline", "/", pmhc_select_out$pmhc)
pmhc_select_out$pmhc <- gsub("colon", ":", pmhc_select_out$pmhc)
pmhc_select_out$pmhc <- gsub("percent", "%", pmhc_select_out$pmhc)
pmhc_select_out$pmhc <- gsub("asterisk", "*", pmhc_select_out$pmhc)
pmhc_select_out$pmhc <- gsub("dash", "-", pmhc_select_out$pmhc)

write.table(pmhc_select_out, paste0("outputs/", args[2], ".pmhc_coverage.tsv"), quote=FALSE, row.names=FALSE, sep='\t')

pmhc_select$pmhc <- gsub("newline", "\n", pmhc_select$pmhc)
pmhc_select$pmhc <- gsub("colon", ":", pmhc_select$pmhc)
pmhc_select$pmhc <- gsub("percent", "%", pmhc_select$pmhc)
pmhc_select$pmhc <- gsub("asterisk", "*", pmhc_select$pmhc)
pmhc_select$pmhc <- gsub("dash", "-", pmhc_select$pmhc)
pmhc_select$pmhc <- gsub("ERV.", "", pmhc_select$pmhc)
print(head(pmhc_select))

p <- ggplot(pmhc_select) + geom_bar(aes(x=reorder(pmhc, cumulative_coverage_percent), y=proportion_cells_expressing, fill=antigen_source), stat="identity") + geom_line(aes(x=reorder(pmhc, cumulative_coverage_percent), y=cumulative_coverage_percent, group=1)) + theme_light() + labs(y="Proportion of Tumor Cells with pMHC", x="pMHCs", title=paste0("Tumor Antigen pMHC Coverage\nAmong ", args[2], " Tumor Cells")) + theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust=0.5))
#    p <- p + guides(fill = "none")
    p <- p + labs(fill="Antigen Source")
    p <- p + theme(legend.title.align=0.5)
    p <- p + theme(axis.text.x = element_text(angle=0, vjust = 0.5, hjust=0.5, size=5))
    p <- p + theme(title=element_text(size=14, face="bold"))
    p <- p + theme(axis.title=element_text(size=12, face="bold"))
#    p <- p + theme(axis.text.x = element_blank())
    ggsave(paste0("plots/", args[2], ".tumor_het.pmhc.pdf"), p)
