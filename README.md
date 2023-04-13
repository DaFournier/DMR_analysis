# DMR_analysis

2019 David Fournier

## Protol for methylation extraction (script not included)
Methylation positions were detected using bismark_methylation_extractor tool with --ignore options set for 3pb in 5' and 5bp in 3', according to fluctuations observed on the methylation extraction profiles at each base position. 

## DMR calling strategy 
To call DMR (differentially methylated regions), we used two different strategies. methylKit, which performs a statistics on the methylation levels, was used with a bin size of 100 and delta change between conditions of 30. For DMR calling with methylKit we used a custom script accessible in the Github repository (https://github.com/tAndreani/DNA-Demethylation-Gadd45/blob/master/Bash/Pipeline.sh). The second strategy is a custom method derived from a previous study of us (Schaefer et al. 2018). This custom method considers a DMR as a region of at least 2 consecutive hemi-CpG with mean change of 30% between WT and condition (either all up or down to detect hyper- or hypoDMR, respectively). Downstream analysis of data and visualization was performed using R and various libraries, including methylKit itself for its collection of useful functions.
