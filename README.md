# AIMsetfinder
Peter Pfaffelhuber, Franziska Grundner-Culemann, Veronika Lipphardt, Franz Baumdicker

Overview:

**AIMsetfinder** is a collection of Rscripts to identify sets of Ancestry Informative Markers (AIMs), that minimize the logloss error of a naive Bayes classifier.

It takes as input:

* SNP data (e.g. 1000 Genomes SNP data or user's own data in **vcf.gz** format)
* biogeographic information (or alternatively any discrete phenotype)

to select a set of specified size of optimal AIMs to classify the samples.

The output is:

* a **vcf.gz** file with the selected AIMs
* a list of SNP identifiers 

which can be used in ancestry inferrence methods.

* Furthermore the posterior probabilites of a naive classifier based on these AIMs are given for the input data.


## Table of contents
  * [Quick start](#quick-start)
  * [Dependencies](#installing-dependencies)
  * [How to run](#how-to-run)
  * [Directory structure and analysis output](#directory-structure-and-analysis-output)
  * [Command line arguments](#command-line-arguments)

### Quick start

```
git clone https://github.com/fbaumdicker/AIMsetfinder.git
cd AIMsetfinder
```

Install [dependencies](#installing-dependencies) and then run the test:
`Rscript example_pipeline.r`

### Installing dependencies

For data analysis as well as for our simulation studies, we rely on
R scripts.
  
#### dependencies
 * For multicore-computing, we require the R-package `parallel`.
 * Since, both data from the 1000 genomes project, which is
  analysed here, and the coalescent simulations, come in vcf-format, we require the
  R-package `vcfR`.
 * For some steps in the analysis, we use `vcftools`
  \cite{Danecek2011} and `bcftools` \cite{Li2011}, both can be installed using
  
```
sudo apt-get install vcftools bcftools
```
#### data resources (optional)

In addition, data from the 1000 genomes project (phase 3) was
downloaded, as well as information on the sampling locations. For
this, we used

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.*
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```
The latter file was renamed {\tt 1000G\_SampleListWithLocations.txt},
and the first row (the header) was removed.


#### for simulations (optional)
  The simulation studies are performed using [msprime](https://github.com/tskit-dev/msprime). Most easily installed via
  `pip3 install msprime`. Msprime is a fast coalescent
  simulator. In particular, structured populations (with varying
  population sizes etc) can be simulated.  We are using the
  python-interface of msprime. See [msprime documentation](https://msprime.readthedocs.io/en/stable/) for more
  information.



#### Overview of dependencies:
  * [msprime](https://github.com/tskit-dev/msprime)
  * [vcftools](https://vcftools.github.io/index.html)
  * [bcftools](https://samtools.github.io/bcftools/bcftools.html)
  * [nnet](https://cran.r-project.org/package=nnet)
  * [vcfR](https://cran.r-project.org/package=vcfR)
  * parallel R package (part of base R)

### How to run
To run the test set: ` Rscript example_pipeline.r `

In `data/sim/ooa/`, you will find the file `ooa_chromosome_1_example.vcf.gz`, a small set of 240 simulated individuals that is used in this tutorial.


### Directory structure and analysis output
The analysis generates the following files:
```
./AIMs                list of identifiers of the chosen AIMs
./AIMs.vcf.gz         corresponding states for all individuals
./predictions.csv     table of posterior probabilites for all classes/BGAs as predicted by naive Bayes using AIMs
./classifications.tab classification into the class with the largest probability as in predictions.csv 
```


In which step different files are produced is described in more details in [readme.pdf](https://github.com/fbaumdicker/AIMsetfinder/blob/master/doc/readme.pdf).

