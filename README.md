# ATACseq Data Cleanup and Analysis
## Overview
- This repository contains the data analysis pipeline for ATAC-seq data comparing the effect of treatment with TBT vs control DMSO across sex and week (week 1 vs week 5). Analyzed changes in chromatin accessibility and found an effect in males from week 1.

## Data Locations on HPC3
- All intermediate files are stored on HPC3 at the following locations:
  - STAR output: `/share/crsp/lab/blumberg/share/STAR_output/trimmed`
  - Sorted BAMSs: `/share/crsp/lab/blumberg/share/sorted_bam/trimmed`
  - Merged BAMs: `/share/crsp/lab/blumberg/share/merged_bams`
  - mtDNA removed BAMs: `/share/crsp/lab/blumberg/share/mtDNAremoved_bams`
  - Deduplicated BAMs: `/share/crsp/lab/blumberg/share/dedup_bam`
  - Broad peaks: `/share/crsp/lab/blumberg/share/broadPeaks` (including bedFiles folder)
  - Narrow peaks: `/share/crsp/lab/blumberg/share/narrowPeaks` (including bedFiles folder)
  - Files generated with R (.RData and .rds): `/share/crsp/lab/blumberg/share/bioinformatics`
  - BigWigs files: `/share/crsp/lab/blumberg/share/bigwigs`

## Steps
* Step 1: `scripts/Step1_trimmed_fastqc_script.sub`
  * Used FASTQC 0.11.9 to generate quality control reports
    * Output is in `trimmed_fastqc` folder
  * Used MultiQC 1.29 to generate report compiling all FASTQC reports
    * Output is `reports/multiqc_report.html`, it is also in `trimmed_fastqc` folder
* Step 2: `scripts/Step2_STAR_samtools_sort_bams.sub`
  * Used STAR 2.7.10a to align FASTQ files with mm10
    * Output is in `/share/crsp/lab/blumberg/share/STAR_output/trimmed` folder in HPC3
  * Used SAMtools 1.15.1 to sort and index BAMs
    * Output is in `/share/crsp/lab/blumberg/share/sorted_bam/trimmed` in HPC3
* Step 3: `scripts/Step3_merge_bams_ATACseq.sub`
  * Used SAMtools 1.15.1 to merge lanes and generate one BAM file per sample * Output is in `/share/crsp/lab/blumberg/share/merged_bams` folder in HPC3
* Step 4: `scripts/Step4_remove_mtDNA.sub`
  * Used SAMtools 1.15.1 to remove mtDNA from BAMs
    * Output is in `/share/crsp/lab/blumberg/share/mtDNAremoved_bams` folder in HPC3
* Step 5: `scripts/Step5_PCR_duplicates_removal.sub`
  * Used PICARD 3.3.0 to remove PCR duplicates from BAMs, and SAMtools 1.15.1 to reindex
    * Output is in `/share/crsp/lab/blumberg/share/dedup_bam` folder in HPC3
* Step 6: `scripts/Step6_MACS2_peakcalling.sub`
  * Used MACS 2.2.7.1 to generate broad and narrow peak files in HPC3, bed files are inside `bedFiles` folder
    * I used the negative control file `MF591-PE.bam`
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
      * Narrow peaks are in `/share/crsp/lab/blumberg/share/narrowPeaks` folder in HPC3, bed files are inside `bedFiles` folder
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
      * Broad peaks are in `/share/crsp/lab/blumberg/share/broadPeaks` folder in HPC3, bed files are inside `bedFiles` folder
* Step 7: `scripts/Step7_flagstat.sub`
  * Used SAMtools 1.15.1 to create flagstat reports
    * Output is `reports/flagstat.txt` and `reports/flagstat_summary.txt`
      * `flagstat.txt` contains flagstat report for each sample
      * `flagstat_summary.txt` contains only the total number of reads in each sample
* Step 8: `bioinformatics/Step8_DiffBind_ATACseq.R`
  * Used R libraries DiffBind and DESeq2 to perform differential analysis
    * R output files are in `/share/crsp/lab/blumberg/share/bioinformatics` in HPC3
      * I ran this on both broad peaks and narrow peaks, everything that was generated with narrow peaks is in `narrow_peaks` folder in HPC3
      * On broad peaks:
        * `diffbind_counted.rds` contains read count at each peak and the consensus peakset
          * 24 Samples, 75642 sites in matrix:
            |  | ID | Condition | Replicate |    Reads | FRiP |
            | :---: | :---: | :---: | :---: | :---: | :---: |
            | 1 | 201fDMSOw1 |  fDMSOw1 |        1 | 1880751 | 0.13 |
            | 2 | 202fDMSOw1 |  fDMSOw1 |        2 | 1014406 | 0.14 |
            | 3 | 203fDMSOw1 |  fDMSOw1 |        3 | 1506817 | 0.12 |
            | 4 | 204fTBTw1  |  fTBTw1 |        1 | 1890111 | 0.11 |
            | 5 | 205fTBTw1  |  fTBTw1 |        2 | 3320502 | 0.11 |
            | 6 | 206fTBTw1  |  fTBTw1 |        3 | 3166823 | 0.10 |
            | 7 | 213mDMSOw1 |  mDMSOw1 |        1 | 3651781 | 0.10 |
            | 8 | 214mDMSOw1 |  mDMSOw1 |        2 | 1010210 | 0.12 |
            | 9 | 215mDMSOw1 |  mDMSOw1 |        3 | 2781559 | 0.10 |
            | 10 | 217mTBTw1 |   mTBTw1 |        1 | 6203853 | 0.12 |
            | 11 | 218mTBTw1 |   mTBTw1 |        2 | 7642028 | 0.11 |
            | 12 | 225fDMSOw5 |  fDMSOw5 |        1 | 4140907 | 0.10 |
            | 13 | 226fDMSOw5 |  fDMSOw5 |        2 | 4052783 | 0.12 |
            | 14 | 227fDMSOw5 |  fDMSOw5 |        3 | 3344052 | 0.13 |
            | 15 | 228fTBTw5 |   fTBTw5 |        1 | 2967667 | 0.14 |
            | 16 | 230fTBTw5 |   fTBTw5 |        2 | 3006112 | 0.11 |
            | 17 | 238mDMSOw5 |  mDMSOw5 |        1 | 11812505 | 0.09 |
            | 18 | 239mDMSOw5 |  mDMSOw5 |        2 | 8544265 | 0.09 |
            | 19 | 240mTBTw5 |   mTBTw5 |        1 | 5687231 | 0.11 |
            | 20 | 241mTBTw5 |   mTBTw5 |        2 | 5310431 | 0.12 |
            | 21 | 242mTBTw5 |   mTBTw5 |        3 | 4835438 | 0.12 |
            | 22 | 297fTBTw1 |   fTBTw1 |        4 | 4047002 | 0.09 |
            | 23 | 298mTBTw1 |   mTBTw1 |        3 | 4103125 | 0.10 |
            | 24 | 299mDMSOw5 |  mDMSOw5 |        3 | 6694580 | 0.07 |
        * `bioinformatics/plots/correlation_affinity.png` is an affinity score-based correlation graph
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/correlation_affinity.png)
        * `diffbind_analyzed.rds` was created after establishing contrasts based on condition (treatment x sex x week) and running DESeq2 analysis with dba.analyze
          * 24 Samples, 73977 sites in matrix:
          |  | ID | Condition | Replicate |    Reads | FRiP |
          | :---: | :---: | :---: | :---: | :---: | :---: |
          | 1 | 201fDMSOw1 |  fDMSOw1 |        1 | 1880751 | 0.12 |
          | 2 | 202fDMSOw1 |  fDMSOw1 |        2 | 1014406 | 0.13 |
          | 3 | 203fDMSOw1 |  fDMSOw1 |        3 | 1506817 | 0.11 |
          | 4 |  204fTBTw1 |   fTBTw1 |        1 | 1890111 | 0.11 |
          | 5 |  205fTBTw1 |   fTBTw1 |        2 | 3320502 | 0.10 |
          | 6 |  206fTBTw1 |   fTBTw1 |        3 | 3166823 | 0.10 |
          | 7 | 213mDMSOw1 |  mDMSOw1 |        1 | 3651781 | 0.10 |
          | 8 | 214mDMSOw1 |  mDMSOw1 |        2 | 1010210 | 0.11 |
          | 9 | 215mDMSOw1 |  mDMSOw1 |        3 | 2781559 | 0.10 |
          | 10 | 217mTBTw1 |   mTBTw1 |        1 | 6203853 | 0.12 |
          | 11 | 218mTBTw1 |   mTBTw1 |        2 | 7642028 | 0.10 |
          | 12 | 225fDMSOw5 |  fDMSOw5 |        1 | 4140907 | 0.09 |
          | 13 | 226fDMSOw5 |  fDMSOw5 |        2 | 4052783 | 0.12 |
          | 14 | 227fDMSOw5 |  fDMSOw5 |        3 | 3344052 | 0.12 |
          | 15 | 228fTBTw5 |   fTBTw5 |        1 | 2967667 | 0.13 |
          | 16 | 230fTBTw5 |   fTBTw5 |        2 | 3006112 | 0.10 |
          | 17 | 238mDMSOw5 |  mDMSOw5 |        1 | 11812505 | 0.08 |
          | 18 | 239mDMSOw5 |  mDMSOw5 |        2 | 8544265 | 0.09 |
          | 19 | 240mTBTw5 |   mTBTw5 |        1 | 5687231 | 0.11 |
          | 20 | 241mTBTw5 |   mTBTw5 |        2 | 5310431 | 0.11 |
          | 21 | 242mTBTw5 |   mTBTw5 |        3 | 4835438 | 0.12 |
          | 22 | 297fTBTw1 |   fTBTw1 |        4 | 4047002 | 0.08 |
          | 23 | 298mTBTw1 |   mTBTw1 |        3 | 4103125 | 0.10 |
          | 24 | 299mDMSOw5 |  mDMSOw5 |        3 | 6694580 | 0.07 |
          
          Design: [~Condition] | 21 Contrasts:
          |  | Factor |  Group | Samples | Group2 | Samples2 | DB.DESeq2 |
          | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
          | 1 | Condition | fDMSOw1 |      3 | fTBTw1 |       4 |        0 |
          | 2 | Condition | fDMSOw1 |      3 | mDMSOw1 |       3 |     1035 |
          | 3 | Condition | fDMSOw1 |      3 | mTBTw1 |       3 |     1021 |
          | 4 | Condition | fDMSOw1 |      3 | fDMSOw5 |       3 |      566 |
          | 5 | Condition | fDMSOw1 |      3 | mDMSOw5 |       3 |     3265 |
          | 6 | Condition | fDMSOw1 |      3 | mTBTw5 |       3 |      662 |
          | 7 | Condition | fTBTw1 |     4 | mDMSOw1 |       3 |      872 |
          | 8 | Condition | fTBTw1 |      4 | mTBTw1 |       3 |        0 |
          | 9 | Condition | fTBTw1 |      4 | fDMSOw5 |       3 |        0 |
          | 10 | Condition | fTBTw1 |      4 | mDMSOw5 |       3 |      789 |
          | 11 | Condition | fTBTw1 |      4 | mTBTw5 |       3 |        7 |
          | 12 | Condition | mDMSOw1 |      3 | mTBTw1 |       3 |     1223 |
          | 13 | Condition | mDMSOw1 |      3 | fDMSOw5 |       3 |     2295 |
          | 14 | Condition | mDMSOw1 |      3 | mDMSOw5 |       3 |     3688 |
          | 15 | Condition | mDMSOw1 |      3 | mTBTw5 |       3 |      376 |
          | 16 | Condition | mTBTw1 |      3 | fDMSOw5 |       3 |        5 |
          | 17 | Condition | mTBTw1 |      3 | mDMSOw5 |       3 |        0 |
          | 18 | Condition | mTBTw1 |      3 | mTBTw5 |       3 |        0 |
          | 19 | Condition | fDMSOw5 |      3 | mDMSOw5 |       3 |      176 |
          | 20 | Condition | fDMSOw5 |      3 | mTBTw5 |       3 |       34 |
          | 21 | Condition | mTBTw5 |      3 | mDMSOw5 |       3 |        0 |
        * `bioinformatics/plots/contrasted_correlations.png` is correlation after applying contrasts
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/contrasted_correlations.png)
        * `bioinformatics/ATAC_DiffPeaks.txt` is a differential binding report on `fTBTw1` vs `mDMSOw1` because this was what the report was on in the HDD
        * `bioinformatics/ATAC_DiffPeaks_malew1.txt` is a differential binding report on `mTBTw1` vs `mDMSOw1`
        * I attempted to create reports comparing `mTBTw5` vs `mDMSOw1`, `fTBTw1` vs `fDMSOw1`, and `fTBTw5` vs `fDMSOw5` but they were not generated
        * `bioinformatics/plots/Venndiagrams_w1.png` is a Venn diagram on all w1 samples
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/Venndiagram_w1.png)
        * `bioinformatics/plots/Venndiagrams_w5.png` is a Venn diagram on all w5 samples
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/Venndiagram_w5.png)
        * `bioinformatics/plots/MAplot_ATAC_malew1.png` is an MA plot on male w1 samples. 1223 peaks had FDR < 0.05
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/MAplot_malew1.png)
        * `bioinformatics/plots/MAplot_ATAC_malew5.png` is an MA plot on male w5 samples. 0 peaks had FDR < 0.05
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/MAplot_ATAC_malew5.png)
        * `bioinformatics/plots/MAplot_ATAC_femalew1.png` is an MA plot on female w1 samples. 0 peaks had FDR < 0.05
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/MAplot_ATAC_femalew1.png)
        * MA plot of female w5 samples was not generated
        * `bioinformatics/plots/PCAplot_ATAC.png` is a PCA plot showing sample clusters by condition
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/PCAplot_ATAC.png)
        * `bioinformatics/plots/boxplot_malew1.png` is a boxplot showing conts per sample a p-values after normalization on male w1 samples
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/boxplot_malew1.png)
          * I attempted to generate on male w5, female w1, and female w5, but they were not generated
        * `bioinformatics/plots/heatmap_ATAC_differential_malew1.png` is a heatmap showing only significant differential peaks
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/heatmap_ATAC_differential_malew1.png)
          * I attempted to generate on male w5, female w1, and female w5, but they were not generated
          * `bioinformatics/plots/heatmap_ATAC_differential_pearson.png` is a pearson correlation heatmap
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/heatmap_ATAC_pearson.png)
          * `bioinformatics/plots/heatmap_ATAC_differential_spearman.png` is a spearman correlation heatmap
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/heatmap_ATAC_spearman.png)
          * I attempted to create a kendall correlation heatmap but the run would crash out in the cluster so I removed that line of code
      * On narrow peaks:
        * `diffbind_counted_narrow.rds` contains read count at each peak and the consensus peakset
          24 Samples, 39860 sites in matrix:
          |  | ID | Condition | Replicate |    Reads | FRiP |
          | :---: | :---: | :---: | :---: | :---: | :---: |
          | 1 | 201fDMSOw1 |  fDMSOw1 |        1 | 1880751 | 0.11 |
          | 2 | 202fDMSOw1 |  fDMSOw1 |        2 | 1014406 | 0.13 |
          | 3 | 203fDMSOw1 |  fDMSOw1 |        3 | 1506817 | 0.11 |
          | 4 |  204fTBTw1 |   fTBTw1 |        1 | 1890111 | 0.11 |
          | 5 |  205fTBTw1 |   fTBTw1 |        2 | 3320502 | 0.09 |
          | 6 |  206fTBTw1 |   fTBTw1 |        3 | 3166823 | 0.09 |
          | 7 | 213mDMSOw1 |  mDMSOw1 |        1 | 3651781 | 0.08 |
          | 8 | 214mDMSOw1 |  mDMSOw1 |        2 | 1010210 | 0.11 |
          | 9 | 215mDMSOw1 |  mDMSOw1 |        3 | 2781559 | 0.09 |
          | 10 | 217mTBTw1 |   mTBTw1 |        1 | 6203853 | 0.10 |
          | 11 | 218mTBTw1 |   mTBTw1 |        2 | 7642028 | 0.09 |
          | 12 | 225fDMSOw5 |   fDMSOw5 |        1 | 4140907 | 0.08 |
          | 13 | 226fDMSOw5 |  fDMSOw5 |        2 | 4052783 | 0.10 |
          | 14 | 227fDMSOw5 |  fDMSOw5 |        3 | 3344052 | 0.11 |
          | 15 | 228fTBTw5 |   fTBTw5 |        1 | 2967667 | 0.11 |
          | 16 | 230fTBTw5 |   fTBTw5 |        2 | 3006112 | 0.09 |
          | 17 | 238mDMSOw5 |  mDMSOw5 |        1 | 11812505 | 0.07 |
          | 18 | 239mDMSOw5 |  mDMSOw5 |        2 | 8544265 | 0.08 |
          | 19 | 240mTBTw5 |   mTBTw5 |        1 | 5687231 | 0.09 |
          | 20 | 241mTBTw5 |   mTBTw5 |        2 | 5310431 | 0.09 |
          | 21 | 242mTBTw5 |   mTBTw5 |        3 | 4835438 | 0.10 |
          | 22 | 297fTBTw1 |   fTBTw1 |        4 | 4047002 | 0.07 |
          | 23 | 298mTBTw1 |   mTBTw1 |        3 | 4103125 | 0.08 |
          | 24 | 299mDMSOw5 |  mDMSOw5 |        3 | 6694580 | 0.05 |
        * `bioinformatics/plots/correlation_affinity_narrow.png` is an affinity score-based correlation graph
![alt text]()
        * `diffbind_analyzed_narrow.rds` was created after establishing contrasts based on condition (treatment x sex x week) and running DESeq2 analysis with dba.analyze
          24 Samples, 38491 sites in matrix:
          |  | ID | Condition | Replicate |    Reads | FRiP |
          | :---: | :---: | :---: | :---: | :---: | :---: |
          | 1 | 201fDMSOw1 |  fDMSOw1 |        1 | 1880751 | 0.10 |
          | 2 | 202fDMSOw1 |  fDMSOw1 |        2 | 1014406 | 0.12 |
          | 3 | 203fDMSOw1 |  fDMSOw1 |        3 | 1506817 | 0.10 |
          | 4 |  204fTBTw1 |   fTBTw1 |        1 | 1890111 | 0.10 |
          | 5 |  205fTBTw1 |   fTBTw1 |        2 | 3320502 | 0.08 |
          | 6 |  206fTBTw1 |   fTBTw1 |        3 | 3166823 | 0.08 |
          | 7 | 213mDMSOw1 |  mDMSOw1 |        1 | 3651781 | 0.08 |
          | 8 | 214mDMSOw1 |  mDMSOw1 |        2 | 1010210 | 0.10 |
          | 9 | 215mDMSOw1 |  mDMSOw1 |        3 | 2781559 | 0.08 |
          | 10 | 217mTBTw1 |   mTBTw1 |        1 | 6203853 | 0.09 |
          | 11 | 218mTBTw1 |   mTBTw1 |        2 | 7642028 | 0.08 |
          | 12 | 225fDMSOw5 |  fDMSOw5 |        1 | 4140907 | 0.08 |
          | 13 | 226fDMSOw5 |  fDMSOw5 |        2 | 4052783 | 0.09 |
          | 14 | 227fDMSOw5 |  fDMSOw5 |        3 | 3344052 | 0.10 |
          | 15 | 228fTBTw5 |   fTBTw5 |        1 | 2967667 | 0.10 |
          | 16 | 230fTBTw5 |   fTBTw5 |        2 | 3006112 | 0.08 |
          | 17 | 238mDMSOw5 |  mDMSOw5 |        1 | 11812505 | 0.06 |
          | 18 | 239mDMSOw5 |  mDMSOw5 |        2 | 8544265 | 0.07 |
          | 19 | 240mTBTw5 |   mTBTw5 |        1 | 5687231 | 0.08 |
          | 20 | 241mTBTw5 |   mTBTw5 |        2 | 5310431 | 0.09 |
          | 21 | 242mTBTw5 |   mTBTw5 |        3 | 4835438 | 0.09 |
          | 22 | 297fTBTw1 |   fTBTw1 |        4 | 4047002 | 0.07 |
          | 23 | 298mTBTw1 |   mTBTw1 |        3 | 4103125 | 0.08 |
          | 24 | 299mDMSOw5 |  mDMSOw5 |        3 | 6694580 | 0.05 |
          
          Design: [~Condition] | 21 Contrasts:
          |  | Factor |  Group | Samples | Group2 | Samples2 | DB.DESeq2 |
          | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
          1 | Condition | fDMSOw1 |      3 | fTBTw1 |       4 |      350 |
          2 | Condition | fDMSOw1 |      3 | mDMSOw1 |       3 |     1118 |
          3 | Condition | fDMSOw1 |      3 | mTBTw1 |       3 |     1014 |
          4 | Condition | fDMSOw1 |      3 | fDMSOw5 |       3 |      511 |
          5 | Condition | fDMSOw1 |      3 | mDMSOw5 |       3 |     3147 |
          6 | Condition | fDMSOw1 |      3 | mTBTw5 |       3 |      696 |
          7 | Condition |  fTBTw1 |      4 | mDMSOw1 |       3 |      796 |
          8 | Condition |  fTBTw1 |      4 | mTBTw1 |       3 |        4 |
          9 | Condition |  fTBTw1 |      4 | fDMSOw5 |       3 |        0 |
          10 | Condition |  fTBTw1 |      4 | mDMSOw5 |       3 |      805 |
          11 | Condition |  fTBTw1 |      4 | mTBTw5 |       3 |       17 |
          12 | Condition | mDMSOw1 |      3 | mTBTw1 |       3 |     1162 |
          13 | Condition | mDMSOw1 |      3 | fDMSOw5 |       3 |     2448 |
          14 | Condition | mDMSOw1 |      3 | mDMSOw5 |       3 |     3503 |
          15 | Condition | mDMSOw1 |      3 | mTBTw5 |       3 |      448 |
          16 | Condition | mTBTw1 |      3 | fDMSOw5 |       3 |        5 |
          17 | Condition | mTBTw1 |      3 | mDMSOw5 |       3 |        0 |
          18 | Condition | mTBTw1 |      3 | mTBTw5 |       3 |        0 |
          19 | Condition | fDMSOw5 |      3 | mDMSOw5 |       3 |      295 |
          20 | Condition | fDMSOw5 |      3 | mTBTw5 |       3 |       39 |
          21 | Condition |  mTBTw5 |      3 | mDMSOw5 |       3 |        0 |
        * `bioinformatics/ATAC_DiffPeaks_narrow.txt` is a differential binding report on `fTBTw1` vs `mDMSOw1` because this was what the report was on in the HDD
        * `bioinformatics/ATAC_DiffPeaks_malew1_narrow.txt` is a differential binding report on `mTBTw1` vs `mDMSOw1`
        * `bioinformatics/ATAC_DiffPeaks_femalew1_narrow.txt` is a differential binding report on `fTBTw1` vs `fDMSOw1`
        * I attempted to create reports comparing `mTBTw5` vs `mDMSOw1`, and `fTBTw5` vs `fDMSOw5` but they were not generated
        * `bioinformatics/plots/Venndiagrams_w1_narrow.png` is a Venn diagram on all w1 samples
        * `bioinformatics/plots/Venndiagrams_w5_narrow.png` is a Venn diagram on all w5 samples
        * `bioinformatics/plots/MAplot_ATAC_malew1_narrow.png` is an MA plot on male w1 samples. 1162 peaks had FDR < 0.05
        * `bioinformatics/plots/MAplot_ATAC_malew5_narrow.png` is an MA plot on male w5 samples. 0 peaks had FDR < 0.05
        * `bioinformatics/plots/MAplot_ATAC_femalew1_narrow.png` is an MA plot on female w1 samples. 350 peaks had FDR < 0.05
        * MA plot of female w5 samples was not generated
        * `bioinformatics/plots/PCAplot_ATAC_narrow.png` is a PCA plot showing sample clusters by condition
        * `bioinformatics/plots/boxplot_malew1_narrow.png` is a boxplot showing conts per sample a p-values after normalization on male w1 samples
        * `bioinformatics/plots/boxplot_femalew1_narrow.png` is a boxplot showing conts per sample a p-values after normalization on female w1 samples
          * I attempted to generate on male w5 and female w5, but they were not generated
        * `bioinformatics/plots/heatmap_ATAC_differential_malew1_narrow.png` is a heatmap showing only significant differential peaks
        * `bioinformatics/plots/heatmap_ATAC_differential_femalew1_narrow.png` is a heatmap showing only significant differential peaks
          * I attempted to generate on male w5 and female w5, but they were not generated
          * `bioinformatics/plots/heatmap_ATAC_differential_pearson_narrow.png` is a pearson correlation heatmap
          * `bioinformatics/plots/heatmap_ATAC_differential_spearman_narrow.png` is a spearman correlation heatmap
          * I attempted to create a kendall correlation heatmap but the run would crash out in the cluster so I removed that line of code
* Step 9: `bioinformatics/Step9_CHIPseeker_ATACseq.R`
  * Used R libraries CHIPseeker and clusterProfiler to perform GO analysis. 
    * R output files are in `bioinformatics` folder. The files generated with narrow peaks are in the `narrow_attempt` folder
      * I uploaded `ATAC_DiffPeaks_malew1.txt` as the peak file and used `annotatePeak` function to match peak coordinates to peak annotations
      * `bioinformatics/plots/piechart_malew1.png` is a pie chart showing genomic distribution
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/piechart_malew1.png)
      * `bioinformatics/plots/upsetplot_malew1.png` shows peak overlaps between annotation categories
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/upsetplot_malew1.png)
      * `bioinformatics/plots/coverageplot_malew1.png` shows tss enrichment
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/coverageplot_malew1.png)
      * `bioinformatics/plots/heatmaptags_ATAC.png` show position relative to TSS
![alt text](https://github.com/glongo001/ATACseq_Batch20/blob/main/bioinformatics/plots/heatmaptags_ATAC.png)
      * `bioinformatics/DiffPeaksGO.txt` is the GO analysis on male w1 samples
