#!/bin/bash

#SBATCH --job-name=HOMERdmr        # job name
##SBATCH --array=1-4
#SBATCH --nodes=1                            # nodes
#SBATCH -p long                             # queue -- worked with long but w. nodeshort benefits from all cores of node
#SBATCH -A imb     #IMB cluster -- imbniehrs on Mogon
#SBATCH -c 10       # nb of cores/task -- 8 worked
#SBTACH --ntasks=1  # 1 task here but helps to allocate more resources
#SBATCH --mem=36000M                          # memory
#SBATCH --time=2-00:00:00                      # time
#SBATCH --mail-user=dafourni@imb-mainz.de  # email
#SBATCH --mail-type=END                      # END: mail when finished; ALL: for everything


##################################
# David Fournier - December 2018 #
##################################

# Objective: find motive enrichement at hyper- and hypoDMR (WT versus Dhe +/-)
# launch with sbatch peaks_HOMER_pipeline.sh


# Modules load

module load bedtools # IMB server version
#module load bio/BEDTools  # Mogon version



# Variables & folders

# adding HOMER folder to the PATH on the fly

PATH=/fsimb/groups/imb-niehrsgr/mallick/dependencies/analysis_tools/homer/bin:/fsimb/groups/imb-niehrsgr/mallick/dependencies/analysis_tools/motif_analysis/gimmemotif/gimmemotifs-0.8.5/src/weblogo:$PATH  # IMB version
#PATH=/project/imbniehrs/mallick/dependencies/homer/bin/:/project/imbniehrs/mallick/dependencies/DownStream_Analysis_Tool/Motif_Analysis_Tool/GimmeMotif/gimmemotifs-0.8.5/src/weblogo:$PATH  # Mogon version

dir="/home/dafourni/folders/Lmna_aging/DMR"
if [ ! -d "${dir}/analysis/HOMER/HOMER_files" ]; then mkdir ${dir}/analysis/HOMER/HOMER_files ; fi #mkdir -p Homer/output/out_${i} ;


# Functions

perform_HOMER_analysis(){  

    fdir=$1  # full path to bed file
    fout=$2 # name to give for output file 
    bg=$3
    
if [ ! -d "${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}" ]; then mkdir ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout} ; fi #mkdir -p Homer/output/out_${i} ;

bedtools getfasta -s -fullHeader -fi /home/dafourni/folders/indexes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed ${dir}/${fdir} -fo ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}.fa 

#${dir}/analysis/HOMER/bash_getfasta.sh ${dir}/${fdir} # get the fasta region for these intervals
#rm ${working_dir}/data_files/data_files/fC_run1summit_vs_peaks_${win}bp_${seq}.bed # keep file is better

if [ -z "$bg" ]; then # if bg is empty
/fsimb/groups/imb-niehrsgr/mallick/dependencies/analysis_tools/homer/bin/scrambleFasta.pl ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}.fa > ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}_scramble.fa # random background
/fsimb/groups/imb-niehrsgr/mallick/dependencies/analysis_tools/homer/bin/findMotifs.pl ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}.fa fasta ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/ -fasta ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}_scramble.fa -p 8 2>  ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/log_${fout}_scramble.txt # IMB cluster version
#/project/imbniehrs/mallick/dependencies/homer/bin/findMotifs.pl ${working_dir}/data_files/fC_run1summit_vs_peaks_${win}bp_${seq}.bed.fa > ${working_dir}/scramble_fC_run1summit_vs_peaks_${win}bp_${seq}.bed.fa${working_dir}/data_files/fC_run1summit_vs_peaks_${win}bp_${seq}.bed.fa ${working_dir}/fasta HOMER_files_5fC_${seq}_${win}pb/ -fasta ${working_dir}/scramble_fC_run1summit_vs_peaks_${win}bp_${seq}.bed.fa -p 8 2> ${working_dir}/log_scramble_fC_run1summit_vs_peaks_${win}bp_${seq}.bed.fa # Mogon version.
else
    #    echo "im in the right place"
#    echo "background file is:" $bg
bedtools getfasta -s -fullHeader -fi /home/dafourni/folders/indexes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed $bg -fo ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}_background_all_hyperDMRs_2CpG.fa
/fsimb/groups/imb-niehrsgr/mallick/dependencies/analysis_tools/homer/bin/findMotifs.pl ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}.fa fasta ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/ -fasta ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}_background_all_hyperDMRs_2CpG.fa -p 8 2>  ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/log_${fout}_background_all_hyperDMRs_2CpG.txt # IMB cluster version
fi

#rm ${dir}/analysis/HOMER/HOMER_files/HOMER_files_${fout}/${fout}_scramble.fa  # TO DO AFTER

##mv ${working_dir}/log_scramble_fC_run1summit_vs_peaks_${win}bp_${seq}.txt  ${working_dir}/HOMER_files/HOMER_files_5fC_${seq}_${win}pb/
#done

}

# MANY HOMER INSTANCES --- MAY REQUIRE TO PARALLELIZE ON MOGON +++

# hyperDMRs
echo "hyperDMRs..."

# methylKit
echo "methylKit..."
echo "bins 100 delta 30..."
perform_HOMER_analysis DMR_methylKit/tables/hyper.bins.100.delta.30.methylKit.bed HyperDMR_methylKit_bins100_delta30 # being redone
echo "bins 100 delta 20..."
perform_HOMER_analysis DMR_methylKit/tables/hyper.bins.100.delta.20.methylKit.bed HyperDMR_methylKit_bins100_delta20 
echo "bins 200 delta 30..."
perform_HOMER_analysis DMR_methylKit/tables/hyper.bins.200.delta.30.methylKit.bed HyperDMR_methylKit_bins200_delta30 
echo "bins 200 delta 20..."
perform_HOMER_analysis DMR_methylKit/tables/hyper.bins.200.delta.20.methylKit.bed HyperDMR_methylKit_bins200_delta20 
echo "bins 200 delta 10..."
perform_HOMER_analysis DMR_methylKit/tables/hyper.bins.200.delta.10.methylKit.bed HyperDMR_methylKit_bins200_delta10 

# MethPipe
echo "MethPipe..."
echo "method1 2CpG..."
perform_HOMER_analysis DMR_MethPipe/method1/HyperDMR.2.CpG.method1.bed HyperDMR_MethPipe1_2CpG
echo "method1 3CpG..."
perform_HOMER_analysis DMR_MethPipe/method1/HyperDMR.3.CpG.method1.bed HyperDMR_MethPipe1_3CpG
echo "method2 2CpG..."
perform_HOMER_analysis DMR_MethPipe/method2/HyperDMR.2.CpG.method2.bed HyperDMR_MethPipe2_2CpG
echo "method2 3CpG..."
perform_HOMER_analysis DMR_MethPipe/method2/HyperDMR.3.CpG.method2.bed HyperDMR_MethPipe2_3CpG
echo "method12 2CpG..."
perform_HOMER_analysis DMR_MethPipe/method_1_2_merged/HyperDMR.2.CpG.method_1_2.bed HyperDMR_MethPipe12_2CpG
echo "methods12 3CpG..."
perform_HOMER_analysis DMR_MethPipe/method_1_2_merged/HyperDMR.3.CpG.method_1_2.bed HyperDMR_MethPipe12_3CpG

# Schaefer et al. 
echo "Schaefer et al..."
echo "2CpG..."
perform_HOMER_analysis DMR_Schaefer_method/hyperDMR.clust2CpG.filtered.bed hyperDMR.Schaefer_2CpG
echo "3CpG..."
perform_HOMER_analysis DMR_Schaefer_method/hyperDMR.clust3CpG.filtered.bed hyperDMR.Schaefer_3CpG


# triple intersection methylKit, MethPipe (two methods merged) and Schaefer

perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hyper.methylKit.MethPipe12.Schaefer.bed hyperDMR.Schaefer_2CpG.MethPipe12_2CpG.methylKit_bin200_delta20
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.methylKit.MethPipe12.Schaefer.bed hyperDMR.Schaefer_3CpG.MethPipe12_2CpG.methylKit_bin200_delta20

# other intersections (2 by 2)

perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hyper.methylKit.MethPipe12.bed hyperDMR.MethPipe12_2CpG.methylKit_bin200_delta20  # same for Schaefer 2 CpG and 3 CpG
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hyper.methylKit.Schaefer.bed hyperDMR.Schaefer_2CpG.methylKit_bin200_delta20
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hyper.MethPipe12.Schaefer.bed hyperDMR.Schaefer_2CpG.MethPipe12_2CpG

perform_HOMER_analysis analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.methylKit.Schaefer.bed hyperDMR.Schaefer_3CpG.methylKit_bin200_delta20
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.MethPipe12.Schaefer.bed hyperDMR.Schaefer_3CpG.MethPipe12_2CpG


# methylKit versus Schaefer et al. for all parameters combinations 

if false;then
for a in 2 3; do 
    for b in 100 200; do 
	     for c in 10 20 30; do
perform_HOMER_analysis analysis/overlaps_Schaefer_methylKit_all_combinations/hyper.Schaefer_${a}CpG.methylKit_bin${b}_delta${c}.bed hyperDMR.Schaefer_${a}CpG.methylKit_bin${b}_delta${c}
	     done
    done
done
fi

# methylKit versus Schaefer et al. for all parameters combinations with p-values associated to Schaefer et al. hyperDMRs

if false;then
a=2
for b in 100 200; do 
    for c in 10 20 30; do
		perform_HOMER_analysis analysis/overlaps_Schaefer_padj_methylKit_all_combinations/hyper.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c}.bed hyperDMR.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c}
		perform_HOMER_analysis analysis/overlaps_Schaefer_padj_methylKit_all_combinations/hyper.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c}.bed hyperDMR.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c} ${dir}/DMR_Schaefer_method/hyperDMR.clust${a}CpG.filtered.bed # use all hyperDMRs as background (3rd argument)
    done
done
fi



# all the peaks found with at least two methods -- TO BE DONE AGAIN -- some regions are included two times

#cat ${dir}/analysis/overlaps_MethPipe12_merged/hyper.methylKit.MethPipe12.bed ${dir}/analysis/overlaps_MethPipe12_merged/hyper.methylKit.Schaefer.bed ${dir}/analysis/overlaps_MethPipe12_merged/hyper.MethPipe12.Schaefer.bed | sortBed -i > ${dir}/analysis/overlaps_MethPipe12_merged/hyper.methylKit.MethPipe12.Schaefer.more.than.one.method.bed
#perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hyper.methylKit.MethPipe12.Schaefer.more.than.one.method.bed hyperDMR.Schaefer_2CpG.MethPipe12_2CpG.methylKit_bin200_delta20.more.than.one.method
#cat ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.methylKit.MethPipe12.bed ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.methylKit.Schaefer.bed ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.MethPipe12.Schaefer.bed | sortBed -i > ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.methylKit.MethPipe12.Schaefer.more.than.one.method.bed
#perform_HOMER_analysis ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hyper.methylKit.MethPipe12.Schaefer.more.than.one.method.bed hyperDMR.Schaefer_3CpG.MethPipe12_2CpG.methylKit_bin200_delta20.more.than.one.method


# hypoDMRs
echo "hypoDMRs..."

# methylKit
echo "methylKit..."
echo "bins 100 delta 30..."
perform_HOMER_analysis DMR_methylKit/tables/hypo.bins.100.delta.30.methylKit.bed HypoDMR_methylKit_bins100_delta30
echo "bins 100 delta 20..."
perform_HOMER_analysis DMR_methylKit/tables/hypo.bins.100.delta.20.methylKit.bed HypoDMR_methylKit_bins100_delta20 
echo "bins 200 delta 30..."
perform_HOMER_analysis DMR_methylKit/tables/hypo.bins.200.delta.30.methylKit.bed HypoDMR_methylKit_bins200_delta30 
echo "bins 200 delta 20..."
perform_HOMER_analysis DMR_methylKit/tables/hypo.bins.200.delta.20.methylKit.bed HypoDMR_methylKit_bins200_delta20 
echo "bins 200 delta 10..."
perform_HOMER_analysis DMR_methylKit/tables/hypo.bins.200.delta.10.methylKit.bed HypoDMR_methylKit_bins200_delta10 

# MethPipe
echo "MethPipe..."
echo "method1 2CpG..."
perform_HOMER_analysis DMR_MethPipe/method1/HypoDMR.2.CpG.method1.bed HypoDMR_MethPipe1_2CpG
echo "method1 3CpG..."
perform_HOMER_analysis DMR_MethPipe/method1/HypoDMR.3.CpG.method1.bed HypoDMR_MethPipe1_3CpG
echo "method2 2CpG..."
perform_HOMER_analysis DMR_MethPipe/method2/HypoDMR.2.CpG.method2.bed HypoDMR_MethPipe2_2CpG
echo "method2 3CpG..."
perform_HOMER_analysis DMR_MethPipe/method2/HypoDMR.3.CpG.method2.bed HypoDMR_MethPipe2_3CpG
echo "method12 2CpG..."
perform_HOMER_analysis DMR_MethPipe/method_1_2_merged/HypoDMR.2.CpG.method_1_2.bed HypoDMR_MethPipe12_2CpG
echo "methods12 3CpG..."
perform_HOMER_analysis DMR_MethPipe/method_1_2_merged/HypoDMR.3.CpG.method_1_2.bed HypoDMR_MethPipe12_3CpG

# Schaefer et al. 
echo "Schaefer et al..."
echo "2CpG..."
perform_HOMER_analysis DMR_Schaefer_method/hypoDMR.clust2CpG.filtered.bed hypoDMR.Schaefer_2CpG
echo "3CpG..."
perform_HOMER_analysis DMR_Schaefer_method/hypoDMR.clust3CpG.filtered.bed hypoDMR.Schaefer_3CpG

# triple intersection methylKit, MethPipe (two methods merged) and Schaefer

perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hypo.methylKit.MethPipe12.Schaefer.bed hypoDMR.Schaefer_2CpG.MethPipe12_2CpG.methylKit_bin200_delta20
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.methylKit.MethPipe12.Schaefer.bed hypoDMR.Schaefer_3CpG.MethPipe12_2CpG.methylKit_bin200_delta20

# other intersections

perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hypo.methylKit.MethPipe12.bed hypoDMR.MethPipe12_2CpG.methylKit_bin200_delta20 # same for Schaefer 2 CpG and 3 CpG
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hypo.methylKit.Schaefer.bed hypoDMR.Schaefer_2CpG.methylKit_bin200_delta20
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged/hypo.MethPipe12.Schaefer.bed hypoDMR.Schaefer_2CpG.MethPipe12_2CpG

perform_HOMER_analysis analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.methylKit.Schaefer.bed hypoDMR.Schaefer_3CpG.methylKit_bin200_delta20
perform_HOMER_analysis analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.MethPipe12.Schaefer.bed hypoDMR.Schaefer_3CpG.MethPipe12_2CpG


# methylKit versus Schaefer et al. for all parameters combinations 

if false;then
for a in 2 3; do 
    for b in 100 200; do 
	     for c in 10 20 30; do
perform_HOMER_analysis analysis/overlaps_Schaefer_methylKit_all_combinations/hypo.Schaefer_${a}CpG.methylKit_bin${b}_delta${c}.bed hypoDMR.Schaefer_${a}CpG.methylKit_bin${b}_delta${c}
	     done
    done
done
fi

# methylKit versus Schaefer et al. for all parameters combinations with p-values associated to Schaefer et al. hypoDMRs

if true;then
a=2
for b in 100 200; do 
    for c in 10 20 30; do
	#	perform_HOMER_analysis analysis/overlaps_Schaefer_padj_methylKit_all_combinations/hypo.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c}.bed hypoDMR.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c}
		perform_HOMER_analysis analysis/overlaps_Schaefer_padj_methylKit_all_combinations/hypo.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c}.bed hypoDMR.Schaefer_${a}CpG_padj.methylKit_bin${b}_delta${c} ${dir}/DMR_Schaefer_method/hypoDMR.clust${a}CpG.filtered.bed # use all hypoDMRs as background (3rd argument)
    done
done
fi


# all the peaks found with at least two methods -- TO BE DONE AGAIN -- some regions are included two times

#cat ${dir}/analysis/overlaps_MethPipe12_merged/hypo.methylKit.MethPipe12.bed ${dir}/analysis/overlaps_MethPipe12_merged/hypo.methylKit.Schaefer.bed ${dir}/analysis/overlaps_MethPipe12_merged/hypo.MethPipe12.Schaefer.bed | sortBed -i > ${dir}/analysis/overlaps_MethPipe12_merged/hypo.methylKit.MethPipe12.Schaefer.more.than.one.method.bed
#perform_HOMER_analysis ${dir}/analysis/overlaps_MethPipe12_merged/hypo.methylKit.MethPipe12.Schaefer.more.than.one.method.bed hypoDMR.Schaefer_2CpG.MethPipe12_2CpG.methylKit_bin200_delta20.more.than.one.method
#cat ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.methylKit.MethPipe12.bed ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.methylKit.Schaefer.bed ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.MethPipe12.Schaefer.bed | sortBed -i > ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.methylKit.MethPipe12.Schaefer.more.than.one.method.bed
#perform_HOMER_analysis ${dir}/analysis/overlaps_MethPipe12_merged_Schaefer3CpG/hypo.methylKit.MethPipe12.Schaefer.more.than.one.method.bed hypoDMR.Schaefer_3CpG.MethPipe12_2CpG.methylKit_bin200_delta20.more.than.one.method
