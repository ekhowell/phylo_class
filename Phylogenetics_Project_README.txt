# NCBI CCDS PRDM9: https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&GO=MainBrowse&DATA=CCDS49963.2
# USCS Genome Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables?position=chr17:15543072-15563331&hgsid=1075670165_3pAqvF1yoQSPd61qGTSNji15a0MT&refGene=pack&refGene_sel=1&refSeqComposite_sel=1&hgFind.matches=NM_144809,
# NCBI PRDM9 Protein Entry: https://www.ncbi.nlm.nih.gov/protein/NP_659058.3
# NCBI PRDM9 mRNA Entry: https://www.ncbi.nlm.nih.gov/nuccore/NM_144809.3
# NCBI Genome Data Viewer PRDM9: https://www.ncbi.nlm.nih.gov/genome/gdv/browser/nucleotide/?id=NM_144809.3
# NCBI PRDM9 Conserved Domains: https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=NP_659058.3

# Notes (TO DO):
# Installing software (vcftools, bcftools...)
# mm10 reference 


# Steps:

#-----------------------------------------
# 1) Obtain VCF file of wild mouse genomes and mm10 reference
#-----------------------------------------

# Go here: http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/

# Download two files (AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz and AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz.tbi)

# For mm10 chr 17 reference, download the PRDM9 FASTA file here: https://genome.ucsc.edu/cgi-bin/hgTables?position=chr17:15543072-15563331&hgsid=1077666803_GdXWz7q9AFAzAPVn9dyr5U7dGFVx&refGene=pack&refGene_sel=1&refSeqComposite_sel=1&hgFind.matches=NM_144809,
# Name file as mm10_prdm9_reference_sequence.fa

#-------------------------------------------
# 2) Obtain genomic coordinates of PRDM9 CDS
#-------------------------------------------

# Get the GTF file
# Use the UCSC Genome Table Browser to obtain a GTF file for PRDM9 in the mm10 reference genome
# The entries into the Table Browser retrieval form should look like this: https://genome.ucsc.edu/cgi-bin/hgTables?position=chr17:15543072-15563331&hgsid=1075670165_3pAqvF1yoQSPd61qGTSNji15a0MT&refGene=pack&refGene_sel=1&refSeqComposite_sel=1&hgFind.matches=NM_144809,
# Download the resulting file as mm10_prdm9_ncbi_refseq_gtf.txt

# Convert the GTF file into a BED file to filter the VCF (only need the CHROM, START, and STOP columns)
cat mm10_prdm9_ncbi_refseq_gtf.txt | grep "NM_144809.2" | grep "CDS" | cut -f1,4,5 > mm10_prdm9_ncbi_refseq_gtf.bed

#---------------------------------------------------
# 3) Extract PRDM9 genomic coordinates from VCF file
#---------------------------------------------------

vcftools --gzvcf AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz --bed mm10_prdm9_ncbi_refseq_gtf.bed --recode --out PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS

#-------------------------------------------------------
# 4) Create FASTA file of PRDM9 sequence for each sample
#-------------------------------------------------------

# Go here to find information about the full VCF file: http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/readme_vcf.txt
# Here is a list of the samples separated out by subspecies and population affiliation: 
Mmd_FRA: Mmd_FRA1_14.variant,Mmd_FRA2_15B.variant,Mmd_FRA3_16B.variant,Mmd_FRA4_18B.variant,Mmd_FRA5_B2C.variant,Mmd_FRA6_C1.variant,Mmd_FRA7_E1.variant,Mmd_FRA8_F1B.variant
Mmd_GER: Mmd_GER1_TP1.variant2,Mmd_GER2_TP121B.variant2,Mmd_GER3_TP17-2.variant2,Mmd_GER4_TP3-92.variant2,Mmd_GER5_TP4a.variant2,Mmd_GER6_TP51D.variant2,Mmd_GER7_TP7-10F1A2.variant2,Mmd_GER8_TP81B.variant2
Mmd_IRA: Mmd_IRA1_AH15.variant3,Mmd_IRA2_AH23.variant3,Mmd_IRA3_JR11.variant3,Mmd_IRA4_JR15.variant3,Mmd_IRA5_JR2-F1C.variant3,Mmd_IRA6_JR5-F1C.variant3,Mmd_IRA7_JR7-F1C.variant3,Mmd_IRA8_JR8-F1A.variant3
Mmd_HEL: Mmd_HEL1_HG06.variant4,Mmd_HEL2_HG08.variant4,Mmd_HEL3_HG13.variant4
Mmm_CZE: Mmm_CZE1_CR12.variant5,Mmm_CZE2_CR13.variant5,Mmm_CZE3_CR16.variant5,Mmm_CZE4_CR23.variant5,Mmm_CZE5_CR25.variant5,Mmm_CZE6_CR29.variant5,Mmm_CZE7_CR46.variant5,Mmm_CZE8_CR59.variant5
Mmm_KAZ: Mmm_KAZ1_AL1.variant6,Mmm_KAZ2_AL16.variant6,Mmm_KAZ3_AL19.variant6,Mmm_KAZ4_AL33.variant6,Mmm_KAZ5_AL38.variant6,Mmm_KAZ6_AL40.variant6,Mmm_KAZ7_AL41.variant6,Mmm_KAZ8_AL42.variant6
Mmm_AFG: Mmm_AFG1_396.variant7,Mmm_AFG2_413.variant7,Mmm_AFG3_416.variant7,Mmm_AFG4_424.variant7,Mmm_AFG5_435.variant7,Mmm_AFG6_444.variant7
Mmc_CAST: Mmc_CAST10_H36.variant8,Mmc_CAST1_H12.variant8,Mmc_CAST2_H14.variant8,Mmc_CAST3_H15.variant8,Mmc_CAST4_H24.variant8,Mmc_CAST5_H26.variant8,Mmc_CAST6_H27.variant8,Mmc_CAST7_H28.variant8,Mmc_CAST8_H30.variant8,Mmc_CAST9_H34.variant8
Ms_SPRE: Ms_SPRE1_SP36.variant9,Ms_SPRE2_SP39.variant9,Ms_SPRE3_SP41.variant9,Ms_SPRE4_SP51.variant9,Ms_SPRE5_SP62.variant9,Ms_SPRE6_SP68.variant9,Ms_SPRE7_SP69.variant9,Ms_SPRE8_SP70.variant9


# BGZIP and TABIX the VCF file
bgzip -c PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf > PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf.gz
tabix -p vcf PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf.gz 

# Run this script to output FASTA files for each sample
bash get_fasta.sh

#---------------------------------------------------------
# 5) Create multiple sequence alignment of PRDM9 sequences
#---------------------------------------------------------

# Fix output fasta

cat Ms_SPRE6_SP68.variant9_prdm9_cds.fa | sed -e '1!{/^>.*/d;}' | grep -v ">" | tr -d "\n\r" | fold -w 50 -s > output_file


# Use something like this to combine the files
for file in *.fa
do
   echo ">$file" >> out.fa
   tail -n +2 $file >> out.fa
   echo >> out.fa
done


# 5) create MSA of FASTA files
# 6) decide on model of evolution to use and maximum likelihood software/parameters to use

# Overview:

# Sequence data aquisition

# Converting to FASTA format

# Performing MSA

# Identifying evolutionary model to use

# Maybe including a partition finder step to distinguish between different PRDM9 domain

# Performing maximum likelihood construction of tree