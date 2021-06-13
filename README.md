This repo contains reproducible workflows for the study *Rapid ecosystem-scale consequences of acute deoxygenation on a Caribbean reef*, also known as Hypocolypse.

Workflows are also available in website form at the following address:

https://hypocolypse.github.io/

Alternatively, you can clone this repo for access to all code.

## Data Availability

For information on all raw data and data products, please see the [Data Availability](https://hypocolypse.github.io/data-availability.html) page. There you will find more details and links to the figshare project site and European Nucleotide Archive project page. Please also check at the bottom of individual workflow pages for access to related data. 

## Workflows

You can access each step of the workflows by using the navigation bar at the top of the [Hypocolypse](https://hypocolypse.github.io/) website or by downloading the R Markdown (.Rmd) files in this repo. Below is a brief description of each workflow, as well as information on how to access the code. Workflows appear in order.

### A. Field Analyses

This section contains information on the various field analyses conducted in the study. Workflows can be found on this [webpage](https://hypocolypse.github.io/field.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/field.Rmd). 

### B. 16S rRNA Analysis

This section contains four separate workflows for processing and analyzing the 16s rRNA data set. 

#### No 1. DADA2 Workflow

In this part we go through the steps of processing raw 16S rRNA read data including assessing read quality, filtering reads, correcting errors, and inferring amplicon sequence variants (ASVs). Workflows can be found on this [webpage](https://hypocolypse.github.io/16s-dada2.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/16s-dada2.Rmd). 

#### No 2. Data Preparation

Next we go through the process of defining sample groups, creating phyloseq objects, removing unwanted samples, and removing contaminant ASVs. Various parts of this section can easily be modified to perform different analyses. For example, if you were only interested in a specific group of samples, you could change the code here to create new phyloseq objects. Workflows can be found on this [webpage](https://hypocolypse.github.io/16s-data-prep.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/16s-data-prep.Rmd). 


#### No 3. Diversity

In this workflow, we compare  the taxonomic diversity of normoxic v. hypoxic samples as well as diversity during and after the event. Workflows can be found on this [webpage](https://hypocolypse.github.io/16s-water.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/16s-water.Rmd). 

#### No 4. Differentially Abundant ASVs

Finally, we wanted to understand how ASVs partitioned between normoxic and hypoxic conditions. We used Indicator Species Analysis (ISA) to identify Differentially Abundant (DA) ASVs across the two oxygen states and then visualized the results. Workflows can be found on this [webpage](https://hypocolypse.github.io/16s-da-asvs.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/16s-da-asvs.Rmd). 

### C. Metagenomics

We sequenced four samples for metagenomic analysis---two from Coral Caye and two from Cayo Roldan---both during and after the event. Coral Caye the control site (not impacted) and Cayo Roldan was the impacted site. This section contains seven workflows. 

#### No 1. Setup Working Environment

For the metagenomic analysis we largely leave the R environment and enter the [anvi'o](http://merenlab.org/software/anvio/) ecosystem. The first major step in our analysis of  metagenomic samples from hypoxic and normoxic conditions is to setup our working environment. We do this in two steps. The first is to install the backbone of our computational operations, which consists of Miniconda and anvio. Next we build the annotation databases and install the associated tools. We provide complete details of these steps. Workflows can be found on this [webpage](https://hypocolypse.github.io/mg-setup.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mg-setup.Rmd). 

#### No 2. Annotation Databases

Next we describe the various database we used to annotate the metagenomes. There are two main types of annotations we are interested in for this metagenomic project---taxonomic and functional---and there are many, many ways to accomplish both of these goals. This next section involves building the databases and installing any additional tools we need for annotation. Workflows can be found on this [webpage](https://hypocolypse.github.io/mg-databases.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mg-databases.Rmd). 

#### No 3. Assembly & Annotations

Here we move to the processing part of the workflow, which involves the following: adapter trimming of raw data; quality filtering of trimmed reads; co-assembling  QCed reads; mapping reads to the assembly; profiling the mapping results; merging profile dbs; classify reads and genes; annotating genes; and running HMM profiles. Workflows can be found on this [webpage](https://hypocolypse.github.io/mg-workflow-1.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mg-workflow-1.Rmd). 

#### No 4. Assembley & Annotation Summary

In this section of the workflow, we assess the QC & assembly results, Kraken short read taxonomy, mapping results, and contig classifications. Workflows can be found on this [webpage](https://hypocolypse.github.io/mg-workflow-2.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mg-workflow-2.Rmd). 


#### No 5. Binning MAGs

In this section of the workflow we reconstruct metagenome assembled genomes (MAGs), first using CONCOCT for automated binning of the assembled contigs followed by manual refinement. Workflows can be found on this [webpage](https://hypocolypse.github.io/mg-binning.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mg-binning.Rmd). 

#### No 6a. MAG02 Phylogenomic Analysis

In this section, we create a workflow to compare MAG02 with publicly available genomes. MAG02 was only found in the hypoxic sample. Workflows can be found on this [webpage](https://hypocolypse.github.io/mag02-phylogenomics.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mag02-phylogenomics.Rmd). 

### No 6b. MAG04 Phylogenomic Analysis

In this section, we create a workflow to compare MAG04 with publicly available genomes. MAG04 was also only found in the hypoxic sample. Workflows can be found on this [webpage](https://hypocolypse.github.io/mag04-phylogenomics.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mag04-phylogenomics.Rmd). 

### No 7. KEGG/KOFam Analysis

The objective of this section is to visualize all of the genes that have KEGG-KOfam annotations and assess their distribution across metagenomes. Workflows can be found on this [webpage](https://hypocolypse.github.io/mg-function.html) or you can download the [Raw Rmd file](https://github.com/hypocolypse/web/blob/master/mg-function.Rmd). 
