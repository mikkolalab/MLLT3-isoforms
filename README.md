# MLLT3-isoforms
Generation and expansion of transplantable human hematopoietic stem cells (HSC) in culture has been hampered by the difficulty of controlling  transitions between different HSC states without losing stem cell identity. Previous studies showed that MLLT3 can promote human HSC expansion in culture.  We hereby document that the function of the full length MLLT3 (MLLT3-L) is modulated by a short isoform (MLLT3-S) that is transcribed from a downstream TSS and lacks the chromatin binding YEATS domain. Unlike MLLT3-L, MLLT3-S did not promote HSC expansion. Nevertheless,  HSCs deficient for either isoform engrafted poorly in immunodeficient mice. Knockdown and overexpression studies revealed opposing effects for MLLT3 isoforms in controlling HSC regulatory programs, including HSC transcription, protein synthesis, oxidative phosphorylation and splicing. While MLLT3-L is expressed during HSC specification, MLLT3-S is induced after HSC emergence and upregulated in fetal liver HSCs, where it suppresses the pro-proliferative fetal gene IGFBP2. Together, MLLT3 isoforms create a stemness operating system that facilitates HSC transitions between expansion and maintenance modes.
All scripts and processed data from this paper could be found at the box directory: https://ucla.box.com/s/9oq50ise3q5nmrlqug68jw8p0xh32m9e

# Bulk RNA sequencing
All scripts and processed data related with Bulk RNA-sequencing data in this paper could be found at the sub-directory "Bulk RNA-sequencing": https://ucla.box.com/s/eud3sr6sqk28qbu4pch1o0bravkeqz98
There are: 1. scripts to generate salmon index (salmon_index.sh) and to perform salmon quantification (salmon_quant.sh); 2. output of salmon quantification (salmon_output); 3. scripts and figures to visulize salmon quantification (quantification); 4. a mapping table to indicate which samples are included in set1, set2, and set4 (set_to_fastq.xlsx).

# Single cell RNA sequencing
The aggreated single cell RNA sequencing Seurat object can be found at the sub-directory "Single-cell RNA sequencing/Seurat object":
https://ucla.box.com/s/6cm7wo6upxy4qlsocs54odw6bw7fuy1s

The scripts used to analyzed the data and generate the figures can be found at the sub-directory "Single-cell RNA sequencing/Seurat object":
https://ucla.box.com/s/9f9fcnlz5ocplakntm5euojo5f47otvq

# PacBio Iso-seq
All scripts and processed data related with PacBio Iso-seq data in this paper could be found at the sub-directory "PacBio Iso-seq": https://ucla.box.com/s/lq3vhki2661kn9d1iol22qrnks8fxqg2

Scripts and processed data related with isoform quantification could be found at "PacBio Iso-seq/Isoform_quantification": https://ucla.box.com/s/zdvife2q56kobnuexdfywfftuzw75aqp. 
There are: 1. raw iso-seq data of one cord blood sample and one fetal liver sample (raw data); 2.scripts to perform adapter removement, polyA and chimeric reads elimination, alignment (s0.raw_read_processing.smk); 3.scripts to perform non-canonical splice junctions correction(s1.process_sams.smk), scripts to perform annotation, novel isoform identification, filteration, and quantification (s2.talon_pipeline.smk, s3.talon_filter.smk); 4. TPM per transript (combined_filt_talon_abundance_filtered.tsv)

Scripts and processed data related with open reading frame could be found at "PacBio Iso-seq/Open_reading_frame_prediction": https://ucla.box.com/s/p0bpykd2jn48lgonuk0dk1bo9bbelccu
There are: 1. Scripts to find open reading frames (subset.sh and find_ORFs.cds.py); 2. Open reading frame predictions with different cutoffs and gene assemblies (custom_gtfs).

Scripts related transcript structure and abundance visulization could be found at "PacBio Iso-seq/Transcripts_plot": : https://ucla.box.com/s/p0bpykd2jn48lgonuk0dk1bo9bbelccu
