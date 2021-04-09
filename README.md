# Phylogenetics Final Report Reproducible Script
#### By: Emma Howell
#### Last Updated: 4/8/21
Purpose: To record steps in phylogenetic analysis of PRDM9 locus in house mouse subspecies complex.

## Note
This README will include all the steps necessary to recreate this analysis. However, all of the required files/data are contained within this repository so there is no need to re-download VCF files or re-fetch GTF and reference files from UCSC's Table Browser.

## Installing Software
Required software: vcftools and bcftools (in a conda environment), MUSCLE, and IQ-TREE
#### Creating Conda Environment
[TO DO]
#### Installing MUSCLE
[TO DO]
#### Installing IQ-TREE
[TO DO]

## Step 1: Get Wild Mouse Genomes
Go here to download the VCF file containing the filtered SNPs of the wild mouse samples (Warning: LARGE)
http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/
Download the files named "AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz" and "AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz.tbi"

## Step 2: Get Genomic Coordinates of PRDM9 Locus
Next, we need a file that tells us what region of the house mouse genome (build mm10) contains the PRDM9 CDS. We can do this using the UCSC Genome Table Browser. This will give us a GTF file for PRDM9 in the mm10 reference genome. The entries into the Table Browser should be filled out like this:
https://genome.ucsc.edu/cgi-bin/hgTables?position=chr17:15543072-15563331&hgsid=1075670165_3pAqvF1yoQSPd61qGTSNji15a0MT&refGene=pack&refGene_sel=1&refSeqComposite_sel=1&hgFind.matches=NM_144809,
Download the output file as "mm10_prdm9_ncbi_refseq_gtf.txt"

The top of the file should look like this:
```
chr17   mm10_ncbiRefSeq stop_codon      15543973        15543975        0.000000        -       .       gene_id "XM_006523988.3"
chr17   mm10_ncbiRefSeq CDS     15543976        15545360        0.000000        -       2       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq exon    15543072        15545360        0.000000        -       .       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq CDS     15553352        15553545        0.000000        -       1       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq exon    15553352        15553545        0.000000        -       .       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq CDS     15554688        15554755        0.000000        -       0       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq exon    15554688        15554755        0.000000        -       .       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq CDS     15555087        15555358        0.000000        -       2       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq exon    15555087        15555358        0.000000        -       .       gene_id "XM_006523988.3"; transc
chr17   mm10_ncbiRefSeq CDS     15555566        15555667        0.000000        -       2       gene_id "XM_006523988.3"; transc
```
We only want the features labeled NM_144809.2 because this is the accession number for the PRDM9 CDS reference sequence. We also need to convert this GTF file into BED file format (only need the CHROM, START, and STOP columns) in order to intersect these genomic coordinates with the records in the VCF file. To do this, we can use the following command:
`cat mm10_prdm9_ncbi_refseq_gtf.txt | grep "NM_144809.2" | grep "CDS" | cut -f1,4,5 > mm10_prdm9_ncbi_refseq_gtf.bed`

Now, you should have a file that looks like this:
```
chr17	15543976	15545360
chr17	15553352	15553545
chr17	15554688	15554755
chr17	15555087	15555358
chr17	15555566	15555667
chr17	15557301	15557457
chr17	15559023	15559072
chr17	15562415	15562522
chr17	15562814	15562937
chr17	15563221	15563301
```

## Step 3: Extract PRDM9 Variants From Whole Genomes
Now that we have a list of the genomic intervals that contain the PRDM9 CDS, we can use this to pull out SNPs in the house mouse genomes that are falling within these intervals.

This vcftools command will output a new VCF file containing only SNPs falling within the regions specified by the BED file:
`vcftools --gzvcf AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz --bed mm10_prdm9_ncbi_refseq_gtf.bed --recode --out PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS`

Go ahead and compress/index this new VCF using the following commands.
```
bgzip -c PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf > PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf.gz
tabix -p vcf PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf.gz 
```

## Step 4: Convert VCF to FASTA
Now we have the variants that overlap our PRDM9 locus, but in order to do a multiple sequence alignment, we need *sequences* in FASTA format, not *variants* in VCF format.

This really just means that we need to "fill in" the space between our variant sites with whatever nucleotide is present in the mm10 reference genome. To do this, we will use the bcftools *consensus* command.

But first, we need a FASTA file that contains our reference sequence for the chunk of genome that is contained within our PRDM9 VCF. To get this reference FASTA, we need to go back to UCSC's Table Browser and this time get the reference sequence for PRDM9's CDS. 

The entries into the Table Browser form should look like this:
https://genome.ucsc.edu/cgi-bin/hgTables?position=chr17:15543072-15563331&hgsid=1077666803_GdXWz7q9AFAzAPVn9dyr5U7dGFVx&refGene=pack&refGene_sel=1&refSeqComposite_sel=1&hgFind.matches=NM_144809,
Make sure the following options are selected in the next window:
https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1077666803_GdXWz7q9AFAzAPVn9dyr5U7dGFVx&hgta_geneSeqType=genomic&hgta_doGenePredSequence=submit
Download the output file as "mm10_prdm9_reference_sequence.fa".

The top of the FASTA file should look something like this:
```
>mm10_ncbiRefSeq_XM_006523988.3_0 range=chr17:15563221-15563301 5'pad=0 3'pad=0 strand=- repeatMasking=none
ATGTCCTGCACCATGAACACCAACAAGCTGGAAGAAAATAGTCCTGAAGA
AGATACAGGGAAATTCGAGTGGAAACCCAAG
>mm10_ncbiRefSeq_XM_006523988.3_1 range=chr17:15562814-15562937 5'pad=0 3'pad=0 strand=- repeatMasking=none
GTCAAAGATGAATTCAAAGACATTTCCATATACTTCTCCAAAGAAGAATG
GGCAGAGATGGGAGAGTGGGAGAAAATTCGCTATAGAAATGTGAAAAGGA
ACTATAAAATGCTGATTTCAATAG
>mm10_ncbiRefSeq_XM_006523988.3_2 range=chr17:15562415-15562522 5'pad=0 3'pad=0 strand=- repeatMasking=none
GTCTCAGAGCCCCTAGACCAGCTTTCATGTGTTACCAAAGGCAAGCAATG
AAACCCCAAATAAATGACAGTGAGGATTCTGATGAAGAGTGGACACCTAA
GCAACAAG
```

Much like with the GTF file, we only want the FASTA records that correspond to the accession number "NM_144809.2". With one record for each exon in the CDS, this should equate to a total of **10** records with the accession number "NM_144809.2" in the header. 

Unfortunately, there's no straightforward way that I've tried yet to automatically extract these specific records [TO DO]. So instead, open up nano and delete the records that do not contain "NM_144809.2" in the header. 

While you're at it, you need to alter the content of the remaining headers so that they are in the format:
`>chr:start-stop`
If the FASTA records for each exon are not in this format in the reference, then bcftools will not be able to recognize them as genomic coordinates when we are converting out VCF to FASTA. Again, there is no better way to do this yet than to just open nano and fix it [TO DO].

(Note: For an edited version of this file, see the refs/ directory). At the end of all this, your new mm10_prdm9_reference_sequence.fa should look like this:
```
>chr17 15543973-15545360
AACTAAGGACAGAAATTCATCCTTGTCTTTTGTGCTCTTTGGCCTTCTCA
AGTCAGAAATTCCTCACTCAACATATGGAATGGAATCATCGCACTGAAAT
CTTCCCAGGAACATCTGCAAGAATAAATCCTAAACCAGGAGATCCCTGTT
CAGATCAGCTTCAGGAACAACATGTTGATTCACAGAACAAAAATGACAAG
GCCAGCAATGAAGTAAAAAGAAAATCCAAACCCAGGCAGAGGATTTCAAC
AACCTTTCCCAGCACACTCAAAGAACAAATGAGATCTGAGGAAAGTAAGA
GAACTGTGGAAGAGCTCAGAACAGGCCAGACAACAAATACAGAGGACACA
GTCAAATCATTTATTGCATCAGAAATCTCAAGTATTGAAAGACAATGTGG
GCAATATTTCAGTGATAAGTCAAATGTCAATGAGCACCAGAAGACACACA
```

Now that we have the reference sequence, we can drop in the variants from our VCF file for each sample in order to create a "consensus" FASTA for each of our sequences. Here is the structure of the bcftools command:
`bcftools consensus -f reference.fa -s sample_name -o output.fa vcf_file.vcf.gz`
Note how we need to specify *which* sample we want to create the FASTA for. This means that this command must be run for *each* sample using the names contained within the VCF file. 

From Harr et al. (2016), here are the names and subpopulation affiliation for each of the sequences contained within the VCF file:
```
Mmd_FRA: Mmd_FRA1_14.variant,Mmd_FRA2_15B.variant,Mmd_FRA3_16B.variant,Mmd_FRA4_18B.variant,Mmd_FRA5_B2C.variant,Mmd_FRA6_C1.variant,Mmd_FRA7_E1.variant,Mmd_FRA8_F1B.variant
Mmd_GER: Mmd_GER1_TP1.variant2,Mmd_GER2_TP121B.variant2,Mmd_GER3_TP17-2.variant2,Mmd_GER4_TP3-92.variant2,Mmd_GER5_TP4a.variant2,Mmd_GER6_TP51D.variant2,Mmd_GER7_TP7-10F1A2.variant2,Mmd_GER8_TP81B.variant2
Mmd_IRA: Mmd_IRA1_AH15.variant3,Mmd_IRA2_AH23.variant3,Mmd_IRA3_JR11.variant3,Mmd_IRA4_JR15.variant3,Mmd_IRA5_JR2-F1C.variant3,Mmd_IRA6_JR5-F1C.variant3,Mmd_IRA7_JR7-F1C.variant3,Mmd_IRA8_JR8-F1A.variant3
Mmd_HEL: Mmd_HEL1_HG06.variant4,Mmd_HEL2_HG08.variant4,Mmd_HEL3_HG13.variant4
Mmm_CZE: Mmm_CZE1_CR12.variant5,Mmm_CZE2_CR13.variant5,Mmm_CZE3_CR16.variant5,Mmm_CZE4_CR23.variant5,Mmm_CZE5_CR25.variant5,Mmm_CZE6_CR29.variant5,Mmm_CZE7_CR46.variant5,Mmm_CZE8_CR59.variant5
Mmm_KAZ: Mmm_KAZ1_AL1.variant6,Mmm_KAZ2_AL16.variant6,Mmm_KAZ3_AL19.variant6,Mmm_KAZ4_AL33.variant6,Mmm_KAZ5_AL38.variant6,Mmm_KAZ6_AL40.variant6,Mmm_KAZ7_AL41.variant6,Mmm_KAZ8_AL42.variant6
Mmm_AFG: Mmm_AFG1_396.variant7,Mmm_AFG2_413.variant7,Mmm_AFG3_416.variant7,Mmm_AFG4_424.variant7,Mmm_AFG5_435.variant7,Mmm_AFG6_444.variant7
Mmc_CAST: Mmc_CAST10_H36.variant8,Mmc_CAST1_H12.variant8,Mmc_CAST2_H14.variant8,Mmc_CAST3_H15.variant8,Mmc_CAST4_H24.variant8,Mmc_CAST5_H26.variant8,Mmc_CAST6_H27.variant8,Mmc_CAST7_H28.variant8,Mmc_CAST8_H30.variant8,Mmc_CAST9_H34.variant8
Ms_SPRE: Ms_SPRE1_SP36.variant9,Ms_SPRE2_SP39.variant9,Ms_SPRE3_SP41.variant9,Ms_SPRE4_SP51.variant9,Ms_SPRE5_SP62.variant9,Ms_SPRE6_SP68.variant9,Ms_SPRE7_SP69.variant9,Ms_SPRE8_SP70.variant9
```

Instead of running the bcftools command an million times, I've thrown those IDs into a list and put the command in a for loop. This "get_fasta.sh" script is contained within the scripts/ directory. (Note: it will need to be run in a directory containing both the VCF and the reference FASTA). Here's how to run it:
`bash get_fasta.sh`

This first bit of the output should look something like this:
```
Applied 0 variants
Applied 0 variants
Applied 0 variants
Applied 0 variants
Applied 3 variants
Applied 3 variants
Applied 0 variants
Applied 3 variants
Applied 0 variants
Applied 0 variants
Applied 0 variants
Applied 0 variants
Applied 0 variants
Applied 3 variants
Applied 0 variants
Applied 0 variants
Applied 4 variants
```
This is bcftools telling us how many sites in the reference FASTA it had to change in order to build each sample FASTA.

## Creating Multi-sample FASTA File
At this point, each sample (of which there are 67) should have its own FASTA file and within each FASTA file, the sequence should be divided up into 10 records that represent the 10 different exons in the PRDM9 CDS.

In order to perform the multiple sequence alignment, we need to merge FASTA records such that there is a single record for the full PRDM9 locus AND we need to pull all of the per-sample FASTA records into one single file to give to MUSCLE. In addition, we want to make sure that we preserve the identity of each sample so we need to include the sample name in the FASTA record. 

To perform this operation, move all of the FASTAs to their own directory (Note: in this repo it is the fastas/ directory) and run the combine_fastas.sh script in that directory. 
`bash combined_fastas.sh`

This should output a single FASTA named "combined_samples.fa". This file should contain 67 record, each of which represents the PRDM9 CDS for a single sample and whose header represents the sample name. Like this:
```
>Mmc_CAST10_H36.variant8_prdm9_cds.fa
AGTCAGAAATTCCTCACTCAACATATGGAATGGAATCATCGCACTGAAAT
CTTCCCAGGAACATCTGCAAGAATAAATCCTAAACCAGGAGATCCCTGTT
CAGATCAGCTTCAGGAACAACATGTTGATTCACAGAACAAAAATGACAAG
GCCAGCAATGAAGTAAAAAGAAAATCCAAACCCAGGCAGAGGATTTCAAC
AACCTTTCCCAGCACACTCAAAGAACAAATGAGATCTGAGGAAAGTAAGA
GAACTGTGGAAGAGCTCAGAACAGGCCAGACAACAAATACAGAGGACACA
GTCAAATCATTTATTGCATCAGAAATCTCAAGTATTGAAAGACAATGTGG
GCAATATTTCAGTGATAAGTCAAATGTCAATGAGCACCAGAAGACACACA
CAGGGGAGAAGCCCTATGTTTGCAGGGAGTGTGGGCGGGGCTTTACACAG
```

Now we are ready for the multiple sequence alignment!