## ----setup, include=FALSE---------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(911)
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(rmarkdown)
library(DT)
library(gridExtra)
library(grid)
knitr::opts_current$get(c(
  "cache",
  "cache.path",
  "cache.rebuild",
  "dependson",
  "autodep"
))

## ----sample_seq_table---------------------
seq_table <- read.table("tables/16s-dada2/sample-seq-info-2019-16s.txt",
                        header = TRUE, sep = "\t")
datatable(seq_table, width = "100%",
          extensions = 'Buttons', options = list(
            scrollX = TRUE,
            dom = 'lfrtipB',
            buttons = c('copy', 'csv', 'excel'),
            pageLength = 5,
            lengthMenu = c(5, 10, 20, 45)
            )
          )

## cutadapt -g {F-ADAPTER} -G {R-ADAPTER} /
##                         -o ${R1}.trimmed.fastq /
##                         -p {R2}.trimmed.fastq /
##                         ${R1} ${R2} /
##                         --discard-untrimmed -e 0.12

## ----list_fastq_filesX------------------
path_X1 <- "01_TRIMMED/RUN01/"
head(list.files(path_X1))


## ----list_fastq_filesY------------------
path_Y2 <- "01_TRIMMED/RUN02/"
head(list.files(path_Y2))

## ----sort_files-------------------------
fnFs_X1 <- sort(list.files(path_X1, pattern = "_R1.trimmed.fastq"))
fnRs_X1 <- sort(list.files(path_X1, pattern = "_R2.trimmed.fastq"))

fnFs_Y2 <- sort(list.files(path_Y2, pattern = "_R1.trimmed.fastq"))
fnRs_Y2 <- sort(list.files(path_Y2, pattern = "_R2.trimmed.fastq"))

## ----split_nameX------------------------
sample.names_X1 <- sapply(strsplit(fnFs_X1, "_"), `[`, 1)
fnFs_X1 <-file.path(path_X1, fnFs_X1)
fnRs_X1 <-file.path(path_X1, fnRs_X1)
sample.names_X1

## ----split_nameY------------------------
sample.names_Y2 <- sapply(strsplit(fnFs_Y2, "_"), `[`, 1)
fnFs_Y2 <-file.path(path_Y2, fnFs_Y2)
fnRs_Y2 <-file.path(path_Y2, fnRs_Y2)
sample.names_Y2


## ----plot_qscoresX---------
p1 <- plotQualityProfile(fnFs_X1[23:30], aggregate = TRUE)
p2 <- plotQualityProfile(fnRs_X1[23:30], aggregate = TRUE)

grid.arrange(p1, p2, nrow = 1)

## ----plot_qscoresY----------
p1 <- plotQualityProfile(fnFs_Y2[1:8], aggregate = TRUE)
p2 <- plotQualityProfile(fnRs_Y2[1:8], aggregate = TRUE)

grid.arrange(p1, p2, nrow = 1)


## ----move_files-------------------------
#Place filtered files in filtered/ subdirectory
filt_path_X1 <- file.path(path_X1, "filtered")
filtFs_X1 <- file.path(filt_path_X1, paste0(sample.names_X1, "_F_filt.fastq.gz"))
filtRs_X1 <- file.path(filt_path_X1, paste0(sample.names_X1, "_R_filt.fastq.gz"))

filt_path_Y2 <- file.path(path_Y2, "filtered")
filtFs_Y2 <- file.path(filt_path_Y2, paste0(sample.names_Y2, "_F_filt.fastq.gz"))
filtRs_Y2 <- file.path(filt_path_Y2, paste0(sample.names_Y2, "_R_filt.fastq.gz"))


## ----filterX----------------------------
out_X1 <- filterAndTrim(fnFs_X1, filtFs_X1, fnRs_X1, filtRs_X1,
                        truncLen=c(260,160), maxN=0, maxEE=c(3,5),
                        truncQ=2, rm.phix=TRUE, compress=TRUE,
                        multithread=TRUE)
head(out_X1)

## ----filterY----------------------------
out_Y2 <- filterAndTrim(fnFs_Y2, filtFs_Y2, fnRs_Y2, filtRs_Y2,
                        truncLen=c(260,160), maxN=0, maxEE=c(5,5),
                        truncQ=2, rm.phix=TRUE, compress=TRUE,
                        multithread=TRUE)
head(out_Y2)

## ----filter_table------------------
datatable(out_X1, width = "80%",
          extensions = 'Buttons', options = list(
            scrollX = TRUE,
            dom = 'lfrtipB',
            buttons = c('copy', 'csv', 'excel'),
            pageLength = 5,
            lengthMenu = c(5, 15, 30, 45)
            )
          )

## ----learn_errors_forwardX---------------
errF_X1 <- learnErrors(filtFs_X1, multithread = TRUE)

## ----learn_errors_forwardY---------------
errF_Y2 <- learnErrors(filtFs_Y2, multithread = TRUE)

## ----learn_errors_reverseX---------------
errR_X1 <- learnErrors(filtRs_X1, multithread = TRUE)

## ----learn_errors_reverseY---------------
errR_Y2 <- learnErrors(filtRs_Y2, multithread = TRUE)

## ----plot_errorFX---------------
plotErrors(errF_X1, nominalQ = TRUE)

## ----plot_errorFY---------------
plotErrors(errF_Y2, nominalQ = TRUE)

## ----plot_errorRX---------------
plotErrors(errR_X1, nominalQ=TRUE)

## ----plot_errorRY---------------
plotErrors(errR_Y2, nominalQ=TRUE)

## ----dereplicate_reads_F---------------
derepFs_X1 <- derepFastq(filtFs_X1)
names(derepFs_X1) <- sample.names_X1

derepFs_Y2 <- derepFastq(filtFs_Y2)
names(derepFs_Y2) <- sample.names_Y2

## ----dereplicate_reads_R---------------
derepRs_X1 <- derepFastq(filtRs_X1)
names(derepRs_X1) <- sample.names_X1

derepRs_Y2 <- derepFastq(filtRs_Y2)
names(derepRs_Y2) <- sample.names_Y2

## ----run_dada2_forwardX1---------------
dadaFs_X1 <- dada(derepFs_X1, err = errF_X1, multithread = TRUE)

## ----run_dada2_forwardY2---------------
dadaFs_Y2 <- dada(derepFs_Y2, err = errF_Y2, multithread = TRUE)

## ----run_dada2_reverseX1---------------
dadaRs_X1 <- dada(derepRs_X1, err = errR_X1, multithread = TRUE)

## ----run_dada2_reverseY2---------------
dadaRs_Y2 <- dada(derepRs_Y2, err = errR_Y2, multithread = TRUE)

## ----inspect_f--------------------------
dadaFs_X1[[1]]

## ----inspect_r--------------------------
dadaRs_X1[[1]]

## ----merge_paired_reads----------------
mergers_X1 <- mergePairs(dadaFs_X1, derepFs_X1, dadaRs_X1, derepRs_X1)
mergers_Y2 <- mergePairs(dadaFs_Y2, derepFs_Y2, dadaRs_Y2, derepRs_Y2)

## ----head_file, eval = FALSE, echo = FALSE----
head(mergers_X1[[1]])

## ----seq_tableX-------------------------
seqtab_X1 <- makeSequenceTable(mergers_X1)
dim(seqtab_X1)

## ----seq_table2X------------------------
table(nchar(getSequences(seqtab_X1)))

## ----seq_tableY-------------------------
seqtab_Y2 <- makeSequenceTable(mergers_Y2)
dim(seqtab_Y2)

## ----seq_table2Y------------------------
table(nchar(getSequences(seqtab_Y2)))

## ----trim_lengthX-----------------------
seqtab_X1.2 <- seqtab_X1[,nchar(colnames(seqtab_X1)) %in% seq(368,383)]
dim(seqtab_X1.2)

## ----trim_length2X---------
table(nchar(getSequences(seqtab_X1.2)))

## ----trim_lengthY-----------------------
seqtab_Y2.2 <- seqtab_Y2[,nchar(colnames(seqtab_Y2)) %in% seq(368,383)]
dim(seqtab_Y2.2)

## ----trim_length2Y---------
table(nchar(getSequences(seqtab_Y2.2)))

## ----saveRDS----------------------------
saveRDS(seqtab_X1.2, "rdata/16s-dada2/seqtab_X1.2.rds")
saveRDS(seqtab_Y2.2, "rdata/16s-dada2/seqtab_Y2.2.rds")

## ----clear_data-----------
remove(list = ls())

## ----read_RDS_files_combo---------------
seqtab.1 <- readRDS("rdata/16s-dada2/seqtab_X1.2.rds")
seqtab.2 <- readRDS("rdata/16s-dada2/seqtab_Y2.2.rds")

st.sum <- mergeSequenceTables(table1 = seqtab.1, table2 = seqtab.2, tables = NULL,
  repeats = "sum", orderBy = "abundance")

## ----save_combo-------------------------
saveRDS(st.sum, "rdata/16s-dada2/combo_run1_run2.rds")

## ----read_combo-------------------------
remove(list = ls())
st.all <- readRDS("rdata/16s-dada2/combo_run1_run2.rds")

## ----chimera_on_ind_runs---------------
seqtab <- removeBimeraDenovo(st.all,
                             method="consensus",
                             multithread=TRUE)
dim(seqtab)

## ----chimera_on_ind_runs2---------------
sum(seqtab)/sum(st.all)

## ----assign_tax_silva, eval = FALSE-----
tax_silva <- assignTaxonomy(
  seqtab, "../taxonomy_databases/silva_nr_v132_train_set.fa.gz",
  multithread = TRUE)

## ----save_image-------------------------
save.image("rdata/16s-dada2/combo_pipeline.rdata")

## ----chimera_on_ind_runsX---------------
#Run01
seqtab_X1.2.nochim <- removeBimeraDenovo(seqtab_X1.2,
                                         method="consensus",
                                         multithread=TRUE)

dim(seqtab_X1.2.nochim)

## ----chimera_on_ind_runsX2---------------
sum(seqtab_X1.2.nochim)/sum(seqtab_X1.2)

## ----chimera_on_ind_runsY---------------
#Run01
seqtab_Y2.2.nochim <- removeBimeraDenovo(seqtab_Y2.2,
                                         method="consensus",
                                         multithread=TRUE)

dim(seqtab_Y2.2.nochim)

## ----chimera_on_ind_runsY2---------------
sum(seqtab_Y2.2.nochim)/sum(seqtab_Y2.2)

## ----build_table_to_track_reads---------------
##Run01
getN_X1 <- function(x) sum(getUniques(x))
track_X1 <-    cbind(out_X1, sapply(dadaFs_X1, getN_X1),
                     sapply(dadaRs_X1, getN_X1), sapply(mergers_X1, getN_X1),
                     rowSums(seqtab_X1.2.nochim))
colnames(track_X1) <- c("input", "filtered", "denoisedF",
                        "denoisedR", "merged", "nonchim")
rownames(track_X1) <- sample.names_X1

##Run02
getN_Y2 <- function(x) sum(getUniques(x))
track_Y2 <-    cbind(out_Y2, sapply(dadaFs_Y2, getN_Y2),
                     sapply(dadaRs_Y2, getN_Y2), sapply(mergers_Y2, getN_Y2),
                     rowSums(seqtab_Y2.2.nochim))
colnames(track_Y2) <- c("input", "filtered", "denoisedF",
                        "denoisedR", "merged", "nonchim")
rownames(track_Y2) <- sample.names_Y2


## ----track_changes----------------------
#Run01
write.table(track_X1, "tables/16s-dada2/RUN01_read_changes.txt",
            sep = "\t", quote = FALSE, col.names=NA)

#Run02
write.table(track_Y2, "tables/16s-dada2/RUN02_read_changes.txt",
            sep = "\t", quote = FALSE, col.names=NA)

