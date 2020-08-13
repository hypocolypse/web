# Working directory

setwd("~/Dropbox (Personal)/STRI-Post-doc/Hypoxia_Bocas/Indicator_species")

# Load SV table
water_seq_table <- read.csv("water_seq_table.csv", stringsAsFactors = FALSE, row.names = 1) # SV table
water_seq_table <- as.data.frame(t(water_seq_table))

# remove all samples except water
# rownames(LEfSe_INPUT_seq_tab)[1:22]->remove
# LEfSe_INPUT_seq_tab <- LEfSe_INPUT_seq_tab[!rownames(LEfSe_INPUT_seq_tab) %in% remove, ]

# Delete columns when sum == 0
water_seq_table <- water_seq_table[, which(colSums(water_seq_table) != 0)]

# write.csv(LEfSe_INPUT_seq_tab, "LEfSe_INPUT_seq_tab_tosee.csv")

# Add First column with name of each group
water_seq_table <- cbind(DO = c("normoxic","normoxic","normoxic","normoxic","normoxic","hypoxic","hypoxic","hypoxic"), water_seq_table)

# install.packages("labdsv")
library(labdsv)

iva <- indval(water_seq_table[,-1], water_seq_table[,1])

gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(water_seq_table[,-1]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]

# Let us see the results
indvalsummary
write.csv(indvalsummary, "indvalsummary.csv")


# Load Tax table
water_tax_table <- read.csv("water_tax_table.csv", stringsAsFactors = FALSE, row.names = 1) # tax table

# Merge indvalsummary with tax names
indvalsummary_tax <- merge(indvalsummary, water_tax_table, by="row.names", all=TRUE)

# Drop rows with NAs
indvalsummary_tax <- indvalsummary_tax[!(is.na(indvalsummary_tax$group)),]

# Rename leves of group factor
# indvalsummary_tax <- as.character(indvalsummary_tax$group)
lapply(indvalsummary_tax, class)
class(indvalsummary_tax$group) = "character"

indvalsummary_tax$group <- ifelse(indvalsummary_tax$group == "1", as.character("hypoxic"), as.character("normoxic"))

write.csv(indvalsummary_tax, "Hypoxia_indvalsummaryWATER_tax.csv", row.names = F)
