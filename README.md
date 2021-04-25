# Phylogenetics Final Report Reproducible Script
#### By: Emma Howell
#### Last Updated: 4/24/21
Purpose: To record the steps used in the phylogenetic analysis of the PRDM9 locus in the house mouse subspecies complex.

## Note
This README will include all the steps necessary to recreate this analysis. However, all of the required input files/data are contained within this repository. This allows users to bypass the steps that involve downloading the (large) VCF files and fetching/formatting the reference files and GTF file obtained from UCSC's Table Browser.

## Installing Software
Required software: vcftools, bcftools, tabix, MUSCLE, and IQ-TREE.

A conda environment will be used to install/manage the packages used in this analysis. Follow the steps below to install Miniconda.

In order to ensure the reproducibility of the analysis, the conda environment used in this pipeline has been packaged into a YAML file contained within the [packages](https://github.com/ekhowell/phylo_class/tree/master/packages) directory. This allows users to bypass the installation steps and ensure that they are using the correct version of each software.

Those wishing to build a conda environment from the YAML file should use the following command.
```
conda env create --file phylo_project_env.yml
```

#### Installing Miniconda
In order to build a conda environment, users must first install Miniconda. Follow the instructions provided [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or see below for a quick guide for MacOS users.

First, download the Linux distribution of Miniconda.
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Then, run the Miniconda installer and follow the prompts.
```
bash Miniconda3-latest-Linux-x86_64.sh
```

#### Creating Conda Environment
Now that Miniconda is installed, we can create a new conda environment called `phylo_project` that will store all of the necessary software.
```
conda create --name phylo_project python=3.8
```

Activate the environment using this command.
```
source activate phylo_project
```

Now we can install the software needed for this project. This includes `vcftools`, `bcftools`, and `tabix` to manipulate the VCF input files. The software needed for the phylogenetic analyses includes `muscle` and `iqtree`.

First, install the software for working with the input files.
```
conda install vcftools
conda install bcftools
conda install tabix
```
Next, install the software for the phylogenetic analyses.
```
conda install muscle
conda install -c bioconda iqtree
```

In order to ensure reproducibility of this analysis (i.e. consistency across software versions) I have packaged this conda environment into the `phylo_project_env.yml` which can be found [here](https://github.com/ekhowell/phylo_class/tree/master/packages).

This command was used to package the built conda environment into a YAML file.
```
conda env export --name phylo_project > phylo_project_env.yml
```

A new conda environment can be built from the YAML file using this command. **NOTE** this step is only neccesary for those that don't want to follow the above steps to install the packages from scratch.
```
conda env create --file phylo_project_env.yml
```

Activate the environment by running the following command.
```
source activate phylo_project
```

## Step 1: Get Wild Mouse Genomes
The sequence data used in this analysis comes from the study published by [Harr et al. (2016)](https://www.nature.com/articles/sdata201675). It contains 67 whole genome sequences sampled collectively from *M. m. domesticus*, *M. m. helgolandicus*, *M. m. musculus*, *M. m. castaneus*, and *M. spretus*. The sequence data is represented in VCF format.

To get the VCF file containing the filtered SNP calls for the wild mouse samples, download the files `AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz` and `AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz.tbi` from the [data repository](http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/). 
```
wget http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz
wget http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz.tbi
```
*WARNING:* These files are **extremely** large. It is recommended that you start with the PRDM9-subset VCF file contained withing the [vcf](https://github.com/ekhowell/phylo_class/tree/master/vcfs) directory.

## Step 2: Get Genomic Coordinates of PRDM9 Locus
Now that we have the whole genome sequence data for the wild mice, we need to extract only the part of the sequence that encapsulates the PRDM9 locus. To do this, we need a file that tells us which region of the house mouse genome (build mm10) contains the coding sequence for PRDM9.

We can do this using the [UCSC Genome Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1095755405_6eTpb3HakfAUV5GkApZw12O9ivSZ) to search for the PRDM9 gene in the mouse genome. This will give us a GTF file containing the location and relevant components of the PRDM9 locus in the mm10 reference genome. The entries into the Table Browser should be filled out as follows:
- Group: `Genes and Gene Predictions`
- Track: `NCBI RefSeq`
- Table: `RefSeq All (ncbiRefSeq)`
- Region: `chr17:15,543,072-15,563,331`

Choose GTF as the file format and download the output table as `mm10_prdm9_ncbi_refseq_gtf.txt`.
**Note:** This GTF file can also be found in the [refs](https://github.com/ekhowell/phylo_class/tree/master/refs) directory.

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

This contains PRDM9 locus information for different tracks in the Table Browser. We only want the features labeled NM_144809.2 because this is the NCBI accession number for the PRDM9 CDS reference sequence. In order to pull these regions out of the VCF file, we also need to convert this GTF file into BED file format. Since the BED file specification only requires a CHR, START, and STOP column, we need to remove the extraneous information contained within the GTF file.

This command will both extract the NM_144809.2 features from the GTF file and remove the additional columns so that the output conforms to the BED file format.
`cat mm10_prdm9_ncbi_refseq_gtf.txt | grep "NM_144809.2" | grep "CDS" | cut -f1,4,5 > mm10_prdm9_ncbi_refseq_gtf.bed`

The resulting file should look like this:
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
**Note:** This final BED file can be found in the [refs](https://github.com/ekhowell/phylo_class/tree/master/refs) directory.

## Step 3: Extract PRDM9 Variants from Whole Genomes
Now that we have the coordinate of the PRDM9 locus CDS contained within the BED file, we can intersect the regions in the BED file with the coordinates in the whole genome VCF file to pull out the SNPs in the house mouse genome that fall within the PRDM9 locus.

This vcftools command will output a new VCF file that contains only those SNPs that overlap with the BED file intervals (i.e. SNPs within the PRDM9 locus):
`vcftools --gzvcf AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz --bed mm10_prdm9_ncbi_refseq_gtf.bed --recode --out PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS`

In order to run the subsequent commands this newly subsetted VCF should be compressed with `bgzip` and indexed with `tabix`.
```
bgzip -c PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf > PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf.gz
tabix -p vcf PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf.gz 
```
**Note:** Both the uncompressed and compressed VCFs can be found within the [vcfs](https://github.com/ekhowell/phylo_class/tree/master/vcfs) directory.

## Step 4: Convert VCF to FASTA
This new VCF file contains the variants that fall within the PRDM9 locus, but in order to perform a multiple sequence alignment and input the data into phylogenetic software, we need *sequences* in FASTA format rather than *variants* in VCF format.

Essentially, this requires "filling in" the space between the variants sites with the non-variants nucleotides that are present in the m10 reference genome. This task can be performed using the `bcftools consensus` command.

In order to run the `consensus` feature, we require a FASTA file that contains the mm10 reference sequence associated with the chunk of genome that we are constructing the consensus for. To obtain the PRDM9 CDS reference sequence, we can once again use the [UCSC Genome Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1095755405_6eTpb3HakfAUV5GkApZw12O9ivSZ).

But first, we need a FASTA file that contains our reference sequence for the chunk of genome that is contained within our PRDM9 VCF. To get this reference FASTA, we need to go back to UCSC's Table Browser and this time get the reference sequence for PRDM9's CDS. 

The entries into the Table Browser should be filled out as follows:
- Group: `Genes and Gene Predictions`
- Track: `NCBI RefSeq`
- Table: `RefSeq All (ncbiRefSeq)`
- Region: `chr17:15,543,079-15,563,323`

Choose sequence as the file format and download the output files as `mm10_prdm9_reference_sequence.fa`.

On the next page, make sure that only the following options are checked in the "Sequence Retrieval Region Options" section:
- CDS Exons
- One FASTA record per region

Check the following options in the "Sequence Formatting Options" sections:
- All upper case
- Mask repeats to lower case

The top of the downloaded FASTA file should look like this:
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

Similar to the GTF file, this FASTA contains sequences for other Table Broswer tracks whereas we only want those corresponding to the accession number NM_144809.2. With one record for each exon in the CDS, there should be a total of **10** records within the downloaded FASTA that contain the accession number NM_144809.2 in the header. 

Unfortunately, there's no straightforward way to automatically extract these specific records. This means that the file must be manually edited. One way to do this is to open the file with `nano` and delete the records that do not contain NM_144809.2 in the header.

In addition, in order to be compatible with the VCF file, the headings of the remaining records beed to be in the format `>chr:start-stop`. Otherwise, `bcftools` will not be able to reconcile the coordinates in the VCF file with the sequences in the FASTA file when building the consensus sequences. These will also require manual editing.

To ensure the reproducibility of this analysis, the manually-edited reference FASTA can be found in the [refs](https://github.com/ekhowell/phylo_class/tree/master/refs) directory. 

The top of the resulting `mm10_prdm9_reference_sequence.fa` file should look like this:
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

Using the reference sequence as a backbone, the variants within the VCF file can be "dropped-in" for each sample in order to create a "filled-in" sequence for each of the mouse samples. Here is the structure of the `bcftools` command used to acheive this:
`bcftools consensus -f reference.fa -s sample_name -o output.fa vcf_file.vcf.gz`
Note how this requires specifying *which* sample the FASTA is being created for. This means that the command must be run for *each* sample using the names contained within the VCF file. 

From [Harr et al. (2016)](https://www.nature.com/articles/sdata201675), below are the names and subpopulation affiliations for each of the sequences contained within the VCF file:
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

Instead of running the `bcftools` command separately for each of the 67 sequences, the `get_fasta.sh` script in the [scripts](https://github.com/ekhowell/phylo_class/tree/master/scripts) directory can be used to automatically run this command for each sample.
`bash get_fasta.sh`
**Note:** This script must be run in a directory containing both the VCF and the reference FASTA file. This may require moving files around.

The top of the screen output should look like this:
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
This represents `bcftools` telling us how many sites in the reference FASTA had to be changed in order to build the FASTA for each sample.

## Creating Multi-sample FASTA File
After running the above script, each of the 67 sequences should have its own FASTA file that represents the PRDM9 locus. Within each FASTA file, the sequence should be divided into 10 records that represent the 10 different exons of the PRDM9 CDS. Some slight formatting needs to be performed before the sequences can be input into MUSCLE for alignment.

In order to perform the multiple sequence alignment, the exon-specific FASTA records within each file must be merged to generate a single record for the entire PRDM9 locus. In addition, the full set of per-sample FASTA sequences must be merged into a single file such that each FASTA record within the combined file represents the full PRDM9 sequence for a given samples. In order to retain information about the subspecies and subpopulation identity of each sequence, the sample name must be included in the header of each record within the combined FASTA file. 

To perform this set of operations, first move all of the FASTAs to their own directory (**Note:** This is the [fastas](https://github.com/ekhowell/phylo_class/tree/master/fastas) directory in the GitHub repository). Next, run the `combine_fastas.sh` script in the directory containing the FASTA files. 
`bash combined_fastas.sh`
**Note:** This script can be found in the [scripts](https://github.com/ekhowell/phylo_class/tree/master/scripts) directory.

This should output a single FASTA file named `combined_samples.fa`. This file should contain 67 records, each of which represents the PRDM9 CDS for a single sample and whose header represents the sample name. The top of the file looks like this:
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
**Note:** This combined FASTA can be found in the [fastas](https://github.com/ekhowell/phylo_class/tree/master/fastas) directory.

## Generating Multiple Sequence Alignment with MUSCLE
To generate the multiple sequence alignment using MUSCLE, run the following command on the FASTA file containing the combined samples.
```
muscle -in combined_samples.fa -phyi -out combined_samples.aln
```
This will produce an alignment in PHYLIP format that is contained within the `combined_samples.aln` file. The top of this output file should look like this:
```
67 2494
Ms_SPRE1_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Ms_SPRE2_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Ms_SPRE3_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Ms_SPRE4_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Ms_SPRE5_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Ms_SPRE6_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Ms_SPRE7_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Ms_SPRE8_S AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
Mmm_AFG2_4 AGTCAGAAAT TCCTCACTCA ACATATGGAA TGGAATCATC GCACTGAAAT
```

This multiple sequence alignment can be used as the starting input for phylogenetic tree estimation using IQ-TREE.

## Generating Phylogeny with IQ-TREE
To perform the phylogenetic tree construction using IQ-TREE, run the following command on the alignment file generated from the MUSCLE step.
```
iqtree -s combined_samples.aln -m MFP -b 100
```

This command strings together multiple of IQ-TREE's functions into a single line of code. The `-m MFP` argument tells IQ-TREE to use ModelFinderPlus to identify the best fit substitution model to use on the data. Then, IQ-TREE constructs a maximum likelihood tree from the alignment supplied by the `-s combined_samples.aln` argument using the best fit substitution model. After constructing the tree, the `-b 100` argument tells IQ-TREE to perform 100 non-parametric bootstrap replicates to evaluate branch supports. 

After running this command, IQ-TREE should direct extensive output to the screen. Once it has finished running, the final section of output should look like this:
```
Analysis results written to: 
  IQ-TREE report:                combined_samples.aln.iqtree
  Maximum-likelihood tree:       combined_samples.aln.treefile
  Likelihood distances:          combined_samples.aln.mldist
  Screen log file:               combined_samples.aln.log

Total CPU time for bootstrap: 172.855 seconds.
Total wall-clock time for bootstrap: 153.781 seconds.

Non-parametric bootstrap results written to:
Bootstrap trees:          combined_samples.aln.boottrees
  Consensus tree:           combined_samples.aln.contree
```
This tells us what is contained in the multiple output files that are produced by IQ-TREE.

The resulting phylogeny contained within the combined_samples.aln.treefile can be viewed and annotated with the [iTOL](https://itol.embl.de/) (interactive Tree of Life) tree viewer and annotation tool.




