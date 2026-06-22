# SRR29831631 Genome Analysis

This repository contains a Galaxy-based genome analysis of SRR29831631, a sequencing dataset labeled as *Clostridium botulinum*. The analysis evaluates whether the genome is consistent with toxigenic *C. botulinum* or more closely related to non-toxigenic *Clostridium sporogenes*.

## Research Question

Is SRR29831631 correctly classified as *Clostridium botulinum*, or is it more genomically consistent with *Clostridium sporogenes*?

## Summary

SRR29831631 was assembled, annotated, screened for botulinum neurotoxin genes, and compared against reference genomes using average nucleotide identity (ANI). The final assembly was high quality, with high completeness and low contamination. No complete `bont` gene was detected, and ANI results showed closer similarity to *C. sporogenes* references than to *C. botulinum* references.

## Main Results

- Final assembly size: 3,984,358 bp
- Scaffolds >=500 bp: 30
- N50: 746,694 bp
- GC content: 27.83%
- CheckM2 completeness: 99.99%
- CheckM2 contamination: 1.66%
- No complete `bont` gene detected
- ANI: >=95% to *C. sporogenes* references and ~92-93% to *C. botulinum* references

## Repository Structure

- `data/final_assembly/`: final filtered scaffold assembly
- `data/raw_accessions/`: SRA and reference genome accession/source notes
- `figures/`: ANI heatmap, PCA, and poster figures
- `methods_notes/`: workflow notes, software versions, and Galaxy dataset key
- `results/`: output tables and reports from QC, assembly assessment, annotation, toxin screening, ANI, and AMR analysis
- `scripts/`: scripts used for figure generation
- `archive_do_not_use/`: older or intermediate files not used in the final analysis

## Workflow Overview

1. Quality control and trimming with FastQC and fastp
2. Genome assembly with SPAdes
3. Assembly filtering to retain scaffolds >=500 bp
4. Assembly quality assessment with QUAST and CheckM2
5. Genome annotation with Prokka and Barrnap
6. Toxin gene screening with ABRicate and BLAST
7. Taxonomic comparison with FastANI
8. Visualization using ANI heatmap and PCA

## Conclusion

Genome-wide ANI comparisons and absence of a complete `bont` gene support interpretation of SRR29831631 as a non-toxigenic *Clostridium sporogenes*-like genome rather than toxigenic *Clostridium botulinum*.

## Data Sources

Raw sequencing data: NCBI SRA SRR29831631  
Reference genomes: NCBI Datasets / NCBI Assembly  
Analysis platform: Galaxy Europe
[Galaxy history: GCC2026](https://usegalaxy.eu/u/priscillagworks/h/gcc2026)

## Author

Priscilla Garcia  
California State University, Stanislaus  
Advisor: Dr. Tricia Van Laar
