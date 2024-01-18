# MLLT3-isoforms
Generation and expansion of transplantable human hematopoietic stem cells (HSC) in culture has been hampered by the difficulty of controlling  transitions between different HSC states without losing stem cell identity. Previous studies showed that MLLT3 can promote human HSC expansion in culture.  We hereby document that the function of the full length MLLT3 (MLLT3-L) is modulated by a short isoform (MLLT3-S) that is transcribed from a downstream TSS and lacks the chromatin binding YEATS domain. Unlike MLLT3-L, MLLT3-S did not promote HSC expansion. Nevertheless,  HSCs deficient for either isoform engrafted poorly in immunodeficient mice. Knockdown and overexpression studies revealed opposing effects for MLLT3 isoforms in controlling HSC regulatory programs, including HSC transcription, protein synthesis, oxidative phosphorylation and splicing. While MLLT3-L is expressed during HSC specification, MLLT3-S is induced after HSC emergence and upregulated in fetal liver HSCs, where it suppresses the pro-proliferative fetal gene IGFBP2. Together, MLLT3 isoforms create a stemness operating system that facilitates HSC transitions between expansion and maintenance modes.
All scripts and processed data from this paper could be found at: https://drive.google.com/drive/folders/1On4yVXCn3LT0wU-WxTyJ9e1rmN0KqBzc?usp=drive_link

# Bulk RNA sequencing
All scripts and processed data related with Bulk RNA-sequencing data in this paper could be found at: https://drive.google.com/drive/u/0/folders/1Ftw7dxFGfhx65MO2PG3GjTKAE2yp2A0Z
There are: 1. scripts to generate salmon index (salmon_index.sh) and tp perform salmon quantification (salmon_quant.sh); 2. output of salmon quantification (salmon_output); 3. scripts and figures to visulize salmon quantification (quantification); 4. a mapping table to indicate which samples are included in set1, set2, and set4 (set_to_fastq.xlsx).

# Single cell RNA sequencing
The aggreated single cell RNA sequencing Seurat object can be found at:
https://drive.google.com/file/d/1_fxYQlQeV_gXRsvVq9YTDhQkPtbUvWuV/view?usp=drive_link

The scripts used to analyzed the data and generate the figures can be found at:
https://drive.google.com/drive/folders/1Qr-px4fJ48p2eh9xx63AKAYJuP_k2jLc

# PacBio Iso-seq
All scripts and processed data related with PacBio Iso-seq data in this paper could be found at: https://drive.google.com/drive/u/0/folders/1ziqMhuJhy0pNcX2cee6zEU-1S9A2HQfT

Scripts and processed data related with isoform quantification could be found at:https://drive.google.com/drive/u/0/folders/1bV1Fro2eKWlT-oYgICXavLCY_W7eF6t3. 
There are: 1. raw iso-seq data of one cord blood sample and one fetal liver sample (raw data); 2.scripts to perform adapter removement, polyA and chimeric reads elimination, alignment (s0.raw_read_processing.smk); 3.scripts to perform non-canonical splice junctions correction(s1.process_sams.smk), scripts to perform annotation, novel isoform identification, filteration, and quantification (s2.talon_pipeline.smk, s3.talon_filter.smk); 4. TPM per transript (combined_filt_talon_abundance_filtered.tsv)

Scripts and processed data related with open reading frame could be found at: https://drive.google.com/drive/u/0/folders/1dheQjZNZN52hKbRRj4tXqDEKvt3CUc70
There are: 1. Scripts to find open reading frames (subset.sh and find_ORFs.cds.py); 2. Open reading frame predictions with different cutoffs and gene assemblies (custom_gtfs).

Scripts related transcript structure and abundance visulization could be found at: https://drive.google.com/drive/u/0/folders/1Fmo5FPaoXO-cBNYSdfDFeGK6O8-GkgoS
