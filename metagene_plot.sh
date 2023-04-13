#!/bin/bash

#SBATCH --job-name=metaPlot 
#SBATCH --nodes=1                            # nb of nodes
#SBATCH -p short                             # queue -- reserves entire node 
#SBATCH -A imb
#SBATCH --cpus-per-task=2         # nb of cores/task
##SBATCH --ntasks=1     # variable set in main.sh script
#SBATCH --mem=24G        # large memory needed for sorting step (8G not enough )
#SBATCH --time=2:00:00                      # time  hh:mm:ss  

# Computes the DNA methylation landscape at a CEBP/hyperDMR site.
# Script by Tommaso Andreani
# adapted by David Fournier January 2019
# original script by TA to be found at the following URL:
# https://raw.githubusercontent.com/tAndreani/DNA-Demethylation-Gadd45/master/Bash/Create.Delta.Mean.Methylation.file.Table.sh
# to be launched using: sbatch metagene_plot.sh


# Module loading - IMB cluster version

module load bedtools
module load datamash
module load deepTools/2.5.1 # version for Debian 8 - IMB cluster runs on Debian 8

# Variables

datadir="/fsimb/groups/imb-niehrsgr/imb_niehrs_2018_DKFZ_XI101_Lmna_WGBS/bismark_methyl_extractions"
outdir="/home/dafourni/folders/Lmna_aging/DMR/analysis/metagene_plots_DMR_CEBP/plots"
if [ ! -d plots ]; then mkdir plots; fi
if [ ! -d intermediary_files ]; then mkdir intermediary_files; fi

###################################################
###Extract methylation values with coverage >= 10##
###then take the common CG shared among all the###
###replicates######################################

if false;then

#From the methylation extracted files (all of them with 43735674 CG extracted ) create a coverage file with methylation values
# +filters CG with Coverage >= 10

for sample in AS-265808-LR-38379 AS-265809-LR-38380 AS-265810-LR-38381 AS-265811-LR-38382; do # sorting step (maybe not needed)
    sort -V -k1,1 -k2,2 ${datadir}/${sample}/${sample}_bismark_bt2_pe.CpG_report.txt > ${sample}_sorted.txt
done
    
for sample in AS-265808-LR-38379 AS-265809-LR-38380 AS-265810-LR-38381 AS-265811-LR-38382; do
    awk '{OFS="\t";if($4+0 > 0 || $5+0 >0 ) print $1"\t"$2-1"\t"$2"\t"$4/($4+$5), $4+$5;}' ${sample}_sorted.txt > tmp.bed    
    cat tmp.bed | awk '$5 >= 10' > ${sample}_cov10.bed 
done


rm tmp.bed

#Create Id values for each CpG
for i in *cov10.bed;
do
    awk -F "\t" '{print $1"_"$2"_"$3}' $i > $i.Id
done 

# sorting step necessary for command "comm" (otherwise comm won't work)
for i in *.Id; do 
sort $i -n > $i.sort
done


#Take common Id

comm -12 AS-265810-LR-38381_cov10.bed.Id.sort  AS-265811-LR-38382_cov10.bed.Id.sort | comm -12 - AS-265808-LR-38379_cov10.bed.Id.sort | comm -12 -  AS-265809-LR-38380_cov10.bed.Id.sort  > Common.CG.all.samples.txt

#Format the Coordinates and sort again
cat Common.CG.all.samples.txt | awk -F "_" '{print $1"\t"$2"\t"$3}' > Common.CG.all.samples.sort.bed
cat Common.CG.all.samples.sort.bed | sort -k 1,1 -k2,2n  > Common.CG.all.samples.sort.as.bed 

#Extract and quality check all are good with the same number
for i in *cov10.bed; do bedtools intersect -a $i -b Common.CG.all.samples.sort.as.bed > $i.common; done


#Create master table
#paste *common | cut -f 1,2,3,4,9,14,19 | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$4"\t"$5}' > Master.Table.wt1.wt2.dhe1.dhe2.txt
paste AS-265808-LR-38379_cov10.bed.common AS-265809-LR-38380_cov10.bed.common AS-265810-LR-38381_cov10.bed.common AS-265811-LR-38382_cov10.bed.common | cut -f 1,2,3,4,9,14,19 > Master.Table.wt1.wt2.dhe1.dhe2.txt 
#paste Sample1_cov10.bed.common Sample2_cov10.bed.common Sample3_cov10.bed.common Sample4_cov10.bed.common | cut -f 1,2,3,4,9,14,19 > Master.Table.wt1.wt2.dhe1.dhe2.txt   # version for TA's pipeline

#Compute Mean Values
cat Master.Table.wt1.wt2.dhe1.dhe2.txt | awk '{print $1"\t"$2"\t"$3"\t"($4+$5)/2"\t"($6+$7)/2}' > Master.Table.wt.dhe.Mean.Values.txt  

#Compute delta
cat Master.Table.wt.dhe.Mean.Values.txt | awk '{print $1"\t"$2"\t"$3"\t"$5-$4}' > Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.txt

#Binning the mouse genome into 100 bp bins
wget https://raw.githubusercontent.com/tAndreani/DNA-Demethylation-Gadd45/master/List.Files.Job.Array/mm10.chr.size
sort -n  mm10.chr.size > tmp.bed
mv tmp.bed mm10.chr.size
bedtools makewindows -g mm10.chr.size -w 100 > mm10.binned.100.bp.ready.bed


##Extract delta mean in bins of 100 bp
sed -i 's/\./,/g' Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.txt # datamash does not support . as floating point
bedtools intersect -a Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.txt -b mm10.binned.100.bp.ready.bed -wb > mm10.binned.100.bp.ready.CpG.mean.values.2.dhe.minus.2.WT

datamash -g 5,6,7 mean 4 < mm10.binned.100.bp.ready.CpG.mean.values.2.dhe.minus.2.WT > Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.txt
sortBed -i Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.txt > Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bedGraph

# bigWig file generation:

#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig -P ~/folders/software/
sed -i 's/,/\./g' Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bedGraph # bedGraphToBigWig does not support , as floating point
~/folders/software/bedGraphToBigWig Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bedGraph mm10.chr.size Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bw

fi

# Cleaning

if false;then
    rm *cov10.bed *.Id *Id.sort Common* *common mm10.binned.100.bp.ready.CpG.mean.values.2.dhe.minus.2.WT Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.txt 
fi    

# metagene plot generation using deepTools


#a="/home/dafourni/folders/Lmna_aging/DMR/analysis/CEBP_ChIP-seq/GSE27826/peaks/SRR125326_CEBPb_sort_peaks_peaks.narrowPeak"
a="/home/dafourni/folders/Lmna_aging/DMR/analysis/CEBP_ChIP-seq/GSE27826/peaks/SRR125326_CEBPb_sort_largepeaks_peaks.broadPeak"
b="/home/dafourni/folders/Lmna_aging/DMR/analysis/overlaps/overlaps_DMR_CEBP/hyper.CEBPb.Schaefer2CpG.bed"
c="/home/dafourni/folders/Lmna_aging/DMR/analysis/overlaps/overlaps_DMR_CEBP/hypo.CEBPb.Schaefer2CpG.bed"
#d="/home/dafourni/folders/Lmna_aging/DMR/analysis/CEBP_ChIP-seq/GSE27826/peaks/SRR125329_CEBPd_sort_peaks_peaks.narrowPeak"
d="/home/dafourni/folders/Lmna_aging/DMR/analysis/CEBP_ChIP-seq/GSE27826/peaks/SRR125329_CEBPd_sort_largepeaks_peaks.broadPeak"
e="/home/dafourni/folders/Lmna_aging/DMR/analysis/overlaps/overlaps_DMR_CEBP/hyper.CEBPd.Schaefer2CpG.bed"
f="/home/dafourni/folders/Lmna_aging/DMR/analysis/overlaps/overlaps_DMR_CEBP/hypo.CEBPd.Schaefer2CpG.bed"

j="Master.Table.wt.dhe.Mean.Values.and.Delta.dhe.minus.WT.bin100bp.bw"

if false;then

#1kb
computeMatrix reference-point --referencePoint center -S ${j} -R ${a} -a 1000 -b 1000 --samplesLabel DNAmet_CEBPBb -out ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join_1kb.gz
computeMatrix reference-point --referencePoint center -S ${j} -R ${b} -a 1000 -b 1000 --samplesLabel DNAmet_CEBPBb -out ${outdir}/matrix_dt_CEBPb_hyperDMR_vs_DNAmet_join_1kb.gz
computeMatrix reference-point --referencePoint center -S ${j} -R ${c} -a 1000 -b 1000 --samplesLabel DNAmet_CEBPBb -out ${outdir}/matrix_dt_CEBPb_hypoDMR_vs_DNAmet_join_1kb.gz
computeMatrix reference-point --referencePoint center -S ${j} -R ${d} -a 1000 -b 1000 --samplesLabel DNAmet_CEBPBd -out ${outdir}/matrix_dt_CEBPd_vs_DNAmet_join_1kb.gz
computeMatrix reference-point --referencePoint center -S ${j} -R ${e} -a 1000 -b 1000 --samplesLabel DNAmet_CEBPBd -out ${outdir}/matrix_dt_CEBPd_hyperDMR_vs_DNAmet_join_1kb.gz
computeMatrix reference-point --referencePoint center -S ${j} -R ${f} -a 1000 -b 1000 --samplesLabel DNAmet_CEBPBd -out ${outdir}/matrix_dt_CEBPd_hypoDMR_vs_DNAmet_join_1kb.gz

#5kb
computeMatrix reference-point --referencePoint center -S ${j} -R ${a} -a 5000 -b 5000 --samplesLabel DNAmet_CEBPBb -out ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join_5kb.gz
#computeMatrix reference-point --referencePoint center -S ${j} -R ${b} -a 5000 -b 5000 --samplesLabel DNAmet_CEBPBb -out ${outdir}/matrix_dt_CEBPb_hyperDMR_vs_DNAmet_join_5kb.gz
#computeMatrix reference-point --referencePoint center -S ${j} -R ${c} -a 5000 -b 5000 --samplesLabel DNAmet_CEBPBb -out ${outdir}/matrix_dt_CEBPb_hypoDMR_vs_DNAmet_join_5kb.gz
#computeMatrix reference-point --referencePoint center -S ${j} -R ${d} -a 5000 -b 5000 --samplesLabel DNAmet_CEBPBd -out ${outdir}/matrix_dt_CEBPd_vs_DNAmet_join_5kb.gz
#computeMatrix reference-point --referencePoint center -S ${j} -R ${e} -a 5000 -b 5000 --samplesLabel DNAmet_CEBPBd -out ${outdir}/matrix_dt_CEBPd_hyperDMR_vs_DNAmet_join_5kb.gz
#computeMatrix reference-point --referencePoint center -S ${j} -R ${f} -a 5000 -b 5000 --samplesLabel DNAmet_CEBPBd -out ${outdir}/matrix_dt_CEBPd_hypoDMR_vs_DNAmet_join_5kb.gz

gunzip ${outdir}/*_1kb.gz
gunzip ${outdir}/*_5kb.gz    
sed -i 's/nan/0/g' ${outdir}/*_1kb
sed -i 's/nan/0/g' ${outdir}/*_5kb
gzip ${outdir}/*1kb
gzip ${outdir}/*5kb

#1kb
plotHeatmap -m ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join_1kb.gz -out ${outdir}/centre_dt_CEBPb_vs_DNAmet_join_1kb.pdf --colorList blue,white,red -x CEBPb_peaks -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
   plotHeatmap -m ${outdir}/matrix_dt_CEBPb_hyperDMR_vs_DNAmet_join_1kb.gz -out ${outdir}/centre_dt_CEBPb_hyperDMR_vs_DNAmet_join_1kb.pdf --colorList blue,white,red -x CEBPb_hyperDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
   plotHeatmap -m ${outdir}/matrix_dt_CEBPb_hypoDMR_vs_DNAmet_join_1kb.gz -out ${outdir}/centre_dt_CEBPb_hypoDMR_vs_DNAmet_join_1kb.pdf --colorList blue,white,red -x CEBPb_hypoDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
   plotHeatmap -m ${outdir}/matrix_dt_CEBPd_vs_DNAmet_join_1kb.gz -out ${outdir}/centre_dt_CEBPd_vs_DNAmet_join_1kb.pdf --colorList blue,white,red -x CEBPd_peaks -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
   plotHeatmap -m ${outdir}/matrix_dt_CEBPd_hyperDMR_vs_DNAmet_join_1kb.gz -out ${outdir}/centre_dt_CEBPd_hyperDMR_vs_DNAmet_join_1kb.pdf --colorList blue,white,red -x CEBPd_hyperDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
   plotHeatmap -m ${outdir}/matrix_dt_CEBPd_hypoDMR_vs_DNAmet_join_1kb.gz -out ${outdir}/centre_dt_CEBPd_hypoDMR_vs_DNAmet_join_1kb.pdf --colorList blue,white,red -x CEBPd_hypoDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no

#5kb
   plotHeatmap -m ${outdir}/matrix_dt_CEBPb_vs_DNAmet_join_5kb.gz -out ${outdir}/centre_dt_CEBPb_vs_DNAmet_join_5kb.pdf --colorList blue,white,red -x CEBPb_peaks -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
#   plotHeatmap -m ${outdir}/matrix_dt_CEBPb_hyperDMR_vs_DNAmet_join_5kb.gz -out ${outdir}/centre_dt_CEBPb_hyperDMR_vs_DNAmet_join_5kb.pdf --colorList blue,white,red -x CEBPb_hyperDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
#   plotHeatmap -m ${outdir}/matrix_dt_CEBPb_hypoDMR_vs_DNAmet_join_5kb.gz -out ${outdir}/centre_dt_CEBPb_hypoDMR_vs_DNAmet_join_5kb.pdf --colorList blue,white,red -x CEBPb_hypoDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
#   plotHeatmap -m ${outdir}/matrix_dt_CEBPd_vs_DNAmet_join_5kb.gz -out ${outdir}/centre_dt_CEBPd_vs_DNAmet_join_5kb.pdf --colorList blue,white,red -x CEBPd_peaks -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
#   plotHeatmap -m ${outdir}/matrix_dt_CEBPd_hyperDMR_vs_DNAmet_join_5kb.gz -out ${outdir}/centre_dt_CEBPd_hyperDMR_vs_DNAmet_join_5kb.pdf --colorList blue,white,red -x CEBPd_hyperDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no
#   plotHeatmap -m ${outdir}/matrix_dt_CEBPd_hypoDMR_vs_DNAmet_join_5kb.gz -out ${outdir}/centre_dt_CEBPd_hypoDMR_vs_DNAmet_join_5kb.pdf --colorList blue,white,red -x CEBPd_hypoDMR -y DNAmet_Delta --refPointLabel center --legendLocation none --boxAroundHeatmaps no

fi

if false;then
   
# plots for CEBPb/d levels at hyper- and hypoDMR

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig -P ~/folders/software/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph -P ~/folders/software/    
~/folders/software/wigToBigWig /home/dafourni/folders/Lmna_aging/DMR/analysis/CEBP_ChIP-seq/GSE27826/CEBPb.wig.gz mm10.chr.size CEBPb.bw
~/folders/software/wigToBigWig /home/dafourni/folders/Lmna_aging/DMR/analysis/CEBP_ChIP-seq/GSE27826/CEBPd.wig.gz mm10.chr.size CEBPd.bw
~/folders/software/bigWigToBedGraph CEBPb.bw  CEBPb.bedGraph # intermediate steps because computeMatrix rejects initial .bw formatting
~/folders/software/bigWigToBedGraph CEBPd.bw  CEBPd.bedGraph
~/folders/software/bedGraphToBigWig CEBPb.bedGraph mm10.chr.size CEBPb.bw
~/folders/software/bedGraphToBigWig CEBPd.bedGraph mm10.chr.size CEBPd.bw

fi

    a="/home/dafourni/folders/Lmna_aging/DMR/DMR_Schaefer_method/hyperDMR.clust2CpG.filtered.bed"
    b="/home/dafourni/folders/Lmna_aging/DMR/DMR_Schaefer_method/hypoDMR.clust2CpG.filtered.bed"
    c="/home/dafourni/folders/Lmna_aging/DMR/analysis/overlaps/overlaps_Schaefer_methylKit_all_combinations/hyper.Schaefer_2CpG.methylKit_bin100_delta30.bed"
    d="/home/dafourni/folders/Lmna_aging/DMR/analysis/overlaps/overlaps_Schaefer_methylKit_all_combinations/hypo.Schaefer_2CpG.methylKit_bin100_delta30.bed"

if false;then
for tf in CEBPb CEBPd; do
	computeMatrix reference-point --referencePoint center -S ${tf}.bw -R ${a} -a 1000 -b 1000 --samplesLabel ${tf}_DNAmet -out ${outdir}/matrix_dt_hyperDMR_Schaefer2CpG_vs_${tf}_1kb.gz
    computeMatrix reference-point --referencePoint center -S ${tf}.bw -R ${b} -a 1000 -b 1000 --samplesLabel ${tf}_DNAmet -out ${outdir}/matrix_dt_hypoDMR_Schaefer2CpG_vs_${tf}_1kb.gz
    computeMatrix reference-point --referencePoint center -S ${tf}.bw -R ${c} -a 1000 -b 1000 --samplesLabel ${tf}_DNAmet -out ${outdir}/matrix_dt_hyperDMR_Schaefer_mKit_vs_${tf}_1kb.gz
    computeMatrix reference-point --referencePoint center -S ${tf}.bw -R ${d} -a 1000 -b 1000 --samplesLabel ${tf}_DNAmet -out ${outdir}/matrix_dt_hypoDMR_Schaefer_mKit_vs_${tf}_1kb.gz
    done
    gunzip ${outdir}/*_1kb.gz    
    sed -i 's/nan/0/g' ${outdir}/*_1kb
    gzip ${outdir}/*1kb
fi

if true;then
# version with 5332 top values (to compare with 5332 hypoDMR:
	computeMatrix reference-point --referencePoint center --sortRegions descend -S CEBPb.bw -R ${a} -a 1000 -b 1000 --samplesLabel CEBPb_DNAmet -out ${outdir}/tmp.gz
	zcat ${outdir}/tmp.gz | head -5332 > tmp
	awk 'NR>1 { print $1 "\t" $2 "\t" $3}' tmp > tmp.bed
	computeMatrix reference-point --referencePoint center --sortRegions descend -S CEBPb.bw -R tmp.bed -a 1000 -b 1000 --samplesLabel CEBPb_DNAmet -out ${outdir}/matrix_dt_hyperDMR_Schaefer2CpG_vs_CEBPb_1kb_top5332.gz
	#rm tmp tmp.bed plots/tmp.gz 

    gunzip ${outdir}/*_top5332.gz    
    sed -i 's/nan/0/g' ${outdir}/*_top5332
    gzip ${outdir}/*_top5332
fi
    
if false;then

  for tf in CEBPb CEBPd; do
    plotHeatmap -m ${outdir}/matrix_dt_hyperDMR_Schaefer2CpG_vs_${tf}_1kb.gz -out ${outdir}/centre_dt_hyperDMR_Schaefer2CpG_vs_${tf}_1kb.pdf --colorList white,red -x HyperDMR -y ${tf}_level --refPointLabel center --legendLocation none --boxAroundHeatmaps no 
    plotHeatmap -m ${outdir}/matrix_dt_hypoDMR_Schaefer2CpG_vs_${tf}_1kb.gz -out ${outdir}/centre_dt_hypoDMR_Schaefer2CpG_vs_${tf}_1kb.pdf --colorList white,red -x HypoDMR -y ${tf}_level --refPointLabel center --legendLocation none --boxAroundHeatmaps no
    plotHeatmap -m ${outdir}/matrix_dt_hyperDMR_Schaefer_mKit_vs_${tf}_1kb.gz -out ${outdir}/centre_dt_hyperDMR_Schaefer_mKit_vs_${tf}_1kb.pdf --colorList white,red -x HyperDMR -y ${tf}_level --refPointLabel center --legendLocation none --boxAroundHeatmaps no
    plotHeatmap -m ${outdir}/matrix_dt_hypoDMR_Schaefer_mKit_vs_${tf}_1kb.gz -out ${outdir}/centre_dt_hypoDMR_Schaefer_mKit_vs_${tf}_1kb.pdf --colorList white,red -x HypoDMR -y ${tf}_level --refPointLabel center --legendLocation none --boxAroundHeatmaps no
    done
fi

if true;then

# version with 5332 top values (to compare with 5332 hypoDMR:	

plotHeatmap -m ${outdir}/matrix_dt_hyperDMR_Schaefer2CpG_vs_CEBPb_1kb_top5332.gz -out ${outdir}/centre_dt_hyperDMR_Schaefer2CpG_vs_CEBPb_1kb_top5332.pdf --colorList white,red -x HypoDMR -y CEBPb_level --refPointLabel center --legendLocation none --boxAroundHeatmaps no

fi

# Put intermediary files in different folder

if false;then
mv *txt intermediary_files/
mv *cov10.bed intermediary_files/
mv *cov10.bed.Id intermediary_files/
mv *Id.sort intermediary_files/
mv *common intermediary_files/
mv Common* intermediary_files/

fi
