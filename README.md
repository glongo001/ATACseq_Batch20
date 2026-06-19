# ATACseq Data Cleanup and Analysis
## Overview
- This repository contains the data analysis pipeline for ATAC-seq data comparing the effect of treatment with TBT vs control DMSO across sex and week (week 1 vs week 5). Analyzed changes in chromatin accessibility and found an effect in males from week 1.

## Data Locations on HPC3
- All intermediate files are stored on HPC3 at the following locations:
  - STAR output: '/share/crsp/lab/blumberg/share/STAR_output/trimmed'
  - Sorted BAMSs: '/share/crsp/lab/blumberg/share/sorted_bam/trimmed'
  - Merged BAMs: '/share/crsp/lab/blumberg/share/merged_bams'
  - mtDNA removed BAMs: '/share/crsp/lab/blumberg/share/mtDNAremoved_bams'
  - Deduplicated BAMs: '/share/crsp/lab/blumberg/share/dedup_bam'
  - Broad peaks: '/share/crsp/lab/blumberg/share/broadPeaks' (including bedFiles folder)
  - Narrow peaks: '/share/crsp/lab/blumberg/share/narrowPeaks' (including bedFiles folder)
  - Files generated with R (.RData and .rds): '/share/crsp/lab/blumberg/share/bioinformatics'
  - BigWigs files: '/share/crsp/lab/blumberg/share/bigwigs'

## Steps
* Step 1: 'scripts/Step1_trimmed_fastqc_script.sub'
  * Used FASTQC 0.11.9 to generate quality control reports
    * Output is in 'trimmed_fastqc' folder
  * Used MultiQC 1.29 to generate report compiling all FASTQC reports
    * Output is 'reports/multiqc_report.html', it is also in 'trimmed_fastqc' folder
* Step 2: 'scripts/Step2_STAR_samtools_sort_bams.sub'
  * Used STAR 2.7.10a to align FASTQ files with mm10
    * Output is in '/share/crsp/lab/blumberg/share/STAR_output/trimmed' folder in HPC3
  * Used SAMtools 1.15.1 to sort and index BAMs
    * Output is in '/share/crsp/lab/blumberg/share/sorted_bam/trimmed' in HPC3
* Step 3: 'scripts/Step3_merge_bams_ATACseq.sub'
  * Used SAMtools 1.15.1 to merge lanes and generate one BAM file per sample * Output is in '/share/crsp/lab/blumberg/share/merged_bams' folder in HPC3
* Step 4: 'scripts/Step4_remove_mtDNA.sub'
  * Used SAMtools 1.15.1 to remove mtDNA from BAMs
    * Output is in '/share/crsp/lab/blumberg/share/mtDNAremoved_bams' folder in HPC3
* Step 5: 'scripts/Step5_PCR_duplicates_removal.sub'
  * Used PICARD 3.3.0 to remove PCR duplicates from BAMs, and SAMtools 1.15.1 to reindex
    * Output is in '/share/crsp/lab/blumberg/share/dedup_bam' folder in HPC3
* Step 6: 'scripts/Step6_MACS2_peakcalling.sub'
  * Used MACS 2.2.7.1 to generate broad and narrow peak files in HPC3, bed files are inside 'bedFiles' folder
    * I used the negative control file 'MF591-PE.bam'
      * I attempted this step without control but it was too strict and lost meaningful data
      * I also attempted it with the downsampled control but it created false peaks, the original negative control file created the cleanest output
    * Parameters for narrow peaks:
      * '''bash
        macs2 callpeak \
          -t "$in_bam" \
          -c "${negative_control_file}" \
          -n "${sample}" \
          -f BAM \
          -g mm \
          -s 75 \
          --outdir "$narrow_path" \
          --nomodel \
          --shift -37 \ #shift reads toward 5' end
          --extsize 73 \ #extend to nucleosome center 73 bp
          -q 0.01 \ #FDR threshold
          --verbose 3
        '''
      * Narrow peaks are in '/share/crsp/lab/blumberg/share/narrowPeaks' folder in HPC3, bed files are inside 'bedFiles' folder
    * Parameters for broad peaks:
      * '''bash
        macs2 callpeak \
          -t "$in_bam" \
          -c "${negative_control_file}" \
          -n "${sample}" \
          -f BAM \
          -g mm \
          -s 75 \
          --outdir "$broad_path" \
          --nomodel \
          --shift -100 \ #shift reads by -100 bp toward 5' end
          --extsize 200 \ #extend to 200 bp
          -q 0.05 \ #FDR threshold, less strict than with narrow peaks
          --broad \ #allow broad peak calling
          --broad-cutoff 0.1 \ #cutoff for broad peak merging
          --verbose 3
        '''
      * Broad peaks are in '/share/crsp/lab/blumberg/share/broadPeaks' folder in HPC3, bed files are inside 'bedFiles' folder
* Step 7: 'scripts/Step7_flagstat.sub'
  * Used SAMtools 1.15.1 to create flagstat reports
    * Output is 'reports/flagstat.txt' and 'reports/flagstat_summary.txt'
      * 'flagstat.txt' contains flagstat report for each sample
      * 'flagstat_summary.txt' contains only the total number of reads in each sample
* Step 8: 'bioinformatics/Step8_DiffBined_ATACseq.R'
  * Used R libraries DiffBind and DESeq2 to perform differential analysis
    * R output files are in '/share/crsp/lab/blumberg/share/bioinformatics' in HPC3
      * I ran this on both broad peaks and narrow peaks, everything that was generated with narrow peaks is in 'narrow_peaks' folder in HPC3
      * On broad peaks:
        * 'diffbind_counted.rds' contains read count at each peak and the consensus peakset
          
        * 'bioinformatics/plots/correlation_affinity.png' is an affinity score-based correlation graph
        * 'diffbind_analyzed.rds' was created after establishing contrasts based on condition (treatment x sex x week) and running DESeq2 analysis with dba.analyze
          
        * 'bioinformatics/ATAC_DiffPeaks.txt' is a differential binding report on 'fTBTw1' vs 'mDMSOw1' because this was what the report was on in the HDD
        * 'bioinformatics/ATAC_DiffPeaks_malew1.txt' is a differential binding report on 'mTBTw1' vs 'mDMSOw1'
        * I attempted to create reports comparing 'mTBTw5' vs 'mDMSOw1', 'fTBTw1' vs 'fDMSOw1', and 'fTBTw5' vs 'fDMSOw5' but they were not generated
        * 'bioinformatics/plots/Venndiagrams_w1.png' is a Venn diagram on all w1 samples
        * 'bioinformatics/plots/Venndiagrams_w5.png' is a Venn diagram on all w5 samples
        * 'bioinformatics/plots/MAplot_ATAC_malew1.png' is an MA plot on male w1 samples. 1223 peaks had FDR < 0.05
        * 'bioinformatics/plots/MAplot_ATAC_malew5.png' is an MA plot on male w5 samples. 0 peaks had FDR < 0.05
        * 'bioinformatics/plots/MAplot_ATAC_femalew1.png' is an MA plot on female w1 samples. 0 peaks had FDR < 0.05
        * MA plot of female w5 samples was not generated
        * 'bioinformatics/plots/PCAplot_ATAC.png' is a PCA plot showing sample clusters by condition
        * 'bioinformatics/plots/boxplot_malew1.png' is a boxplot showing conts per sample a p-values after normalization on male w1 samples
          * I attempted to generate on male w5, female w1, and female w5, but they were not generated
        * 'bioinformatics/plots/heatmap_ATAC_differential_malew1.png' is a heatmap showing only significant differential peaks
          * I attempted to generate on male w5, female w1, and female w5, but they were not generated
          * 'bioinformatics/plots/heatmap_ATAC_differential_pearson.png' is a pearson correlation heatmap
          * 'bioinformatics/plots/heatmap_ATAC_differential_spearman.png' is a spearman correlation heatmap
          * I attempted to create a kendall correlation heatmap but the run would crash out in the cluster so I removed that line of code
      * On narrow peaks:
        * 'diffbind_counted_narrow.rds' contains read count at each peak and the consensus peakset
          
        * 'bioinformatics/plots/correlation_affinity_narrow.png' is an affinity score-based correlation graph
        * 'diffbind_analyzed_narrow.rds' was created after establishing contrasts based on condition (treatment x sex x week) and running DESeq2 analysis with dba.analyze
          
        * 'bioinformatics/ATAC_DiffPeaks_narrow.txt' is a differential binding report on 'fTBTw1' vs 'mDMSOw1' because this was what the report was on in the HDD
        * 'bioinformatics/ATAC_DiffPeaks_malew1_narrow.txt' is a differential binding report on 'mTBTw1' vs 'mDMSOw1'
        * I attempted to create reports comparing 'mTBTw5' vs 'mDMSOw1', 'fTBTw1' vs 'fDMSOw1', and 'fTBTw5' vs 'fDMSOw5' but they were not generated
        * 'bioinformatics/plots/Venndiagrams_w1_narrow.png' is a Venn diagram on all w1 samples
        * 'bioinformatics/plots/Venndiagrams_w5_narrow.png' is a Venn diagram on all w5 samples
        * 'bioinformatics/plots/MAplot_ATAC_malew1_narrow.png' is an MA plot on male w1 samples. 1223 peaks had FDR < 0.05
        * 'bioinformatics/plots/MAplot_ATAC_malew5_narrow.png' is an MA plot on male w5 samples. 0 peaks had FDR < 0.05
        * 'bioinformatics/plots/MAplot_ATAC_femalew1_narrow.png' is an MA plot on female w1 samples. 0 peaks had FDR < 0.05
        * MA plot of female w5 samples was not generated
        * 'bioinformatics/plots/PCAplot_ATAC_narrow.png' is a PCA plot showing sample clusters by condition
        * 'bioinformatics/plots/boxplot_malew1_narrow.png' is a boxplot showing conts per sample a p-values after normalization on male w1 samples
          * I attempted to generate on male w5, female w1, and female w5, but they were not generated
        * 'bioinformatics/plots/heatmap_ATAC_differential_malew1_narrow.png' is a heatmap showing only significant differential peaks
          * I attempted to generate on male w5, female w1, and female w5, but they were not generated
          * 'bioinformatics/plots/heatmap_ATAC_differential_pearson_narrow.png' is a pearson correlation heatmap
          * 'bioinformatics/plots/heatmap_ATAC_differential_spearman_narrow.png' is a spearman correlation heatmap
          * I attempted to create a kendall correlation heatmap but the run would crash out in the cluster so I removed that line of code
        * 
