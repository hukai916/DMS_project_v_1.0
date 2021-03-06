# Intro
This pipeline integrates several bioinformatics tools/packages with the aim of facilitating the analysis of variant enrichment analysis. Two approaches (command line and GUI) are provided. <br>
The statistical enrichment score is derived from Enrich2 results. 6 enrichment plots will be generated by default including 3 amino acids based and 3 codon based figures. Details regarding these plots can be found in Result session below.
This pipeline is tested on Mac OS X (10.13.6) with 16G memory.

# Usage
## Setup
We first need to install the required packages and configure the environment by running the setup.sh. <br>
Download and unzip the git folder, open terminal, then: <br>
**cd scripts/** <br>
**bash setup.sh** <br>
It will take several minutes to install all components.

## GUI approach
After setup, do: <br>
**bash gui.sh** <br>
Input wild type sequence, raw reads file and mutation sites information accordingly.

## Command line approach
To use command line approach, we need to provide with a configuration file. Then do: <br>
**bash wrapper.sh configuration.txt** <br>
The configuration_template.txt file can be found inside script folder.

# Result
The results will be saved in plot_output folder. <br>
By default, 6 plots will be created: 3 are amino acid based and 3 are codon based. <br>
The first plot is raw figure containing all amino acid/codon columns and standard error, the second and third are simplified version of raw figure with empty columns removed and/or standard error bar removed. <br>
Sample output plots are provided in sample_plot_output folder.

**Question and bug report to Kai Hu: (kxh365@@@psu.edu).** <br>
