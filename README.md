# Investigating genome size variation with k-mers

A repository accompanying the manuscript ["Measuring the invisible – The sequences causal of genome size differences in eyebrights (*Euphrasia*) revealed by k-mers"](https://doi.org/10.1101/2021.11.09.467866).

## Motivation
Intraspecific genome size (GS) variation may be due to "presence/absence variation" in single-copy regions or in genomic repeats. To-date, studies targetting the sequence underpinning GS variation commonly use low pass sequencing ("genome skimming") data analysed with the [RepeateExplorer](https://repeatexplorer-elixir.cerit-sc.cz) pipeline. Such studies have convincingly identified repeats involved in GS variation, but they necessarily paint an incomplete picture. Using genome skimming data, it is not possible to assess the contribution to GS variation of low and single-copy sequences.

Here, we use an alternative approach. We use k-mers (short sub-sequences of length 21 generated from sequencing reads) from high-coverage sequencing data sets. We compare k-mer inventories between individuals, which allows us the assess the role of all genomic copy-number classes, from single-copy sequences to highly repetitive satellite DNAs.

## Requirements and dependencies
1. **K-mer tool kit.** You will need to have [KMC3](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) installed. KMC3 can be set up with anaconda, for instance by running `conda install -c bioconda kmc`, or (generating a new environment) `conda create -n kmc -c bioconda kmc`.
2. **File links.** You should rename (or generate links to) your sequencing data files so that each sample has a unique prefix that can be used to easily select all of an individual's files.
3. **Quality filtering/trimming.** (Optionally, but recommended) trim and clean your sequencing data. Sequnecing errors do not matter much. They generate unique k-mers that do not significantly affect estimates. Howerver, sequencing adapter contaminations can show as high-copy number k-mers, biasing genome size estimates. We used [fastp](https://github.com/OpenGene/fastp).
4. **Oraganellar assemblies.** If possible, get sequences for the plastid and mitochrondrial genome. You may choose to assemble *de novo* from your data using [GetOrganelle](https://github.com/Kinggerm/GetOrganelle). These assemblies can then be used to remove organelle k-mers from your data, which would otherwise bias genome size estimates. 

## Running the pipeline
The pipeline has three steps:
1. Generation of k-mer databases and k-mer spectra. These need to be analysed manually (for instance using [Tetmer](https://github.com/hannesbecher/shiny-k-mers)) to assess the sequnecing coverage.
2. Generation of the scaled and binned joint k-mer spectra

### In more detail


