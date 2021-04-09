#!/bin/bash
# Purpose: To create consensus sequence FASTA for each sample

# Create list of sample names as specified in VCF file
declare -a StringArray=("Mmd_FRA1_14.variant" "Mmd_FRA2_15B.variant" "Mmd_FRA3_16B.variant" "Mmd_FRA4_18B.variant" "Mmd_FRA5_B2C.variant" "Mmd_FRA6_C1.variant" "Mmd_FRA7_E1.variant" "Mmd_FRA8_F1B.variant" "Mmd_GER1_TP1.variant2" "Mmd_GER2_TP121B.variant2" "Mmd_GER3_TP17-2.variant2" "Mmd_GER4_TP3-92.variant2" "Mmd_GER5_TP4a.variant2" "Mmd_GER6_TP51D.variant2" "Mmd_GER7_TP7-10F1A2.variant2" "Mmd_GER8_TP81B.variant2" "Mmd_IRA1_AH15.variant3" "Mmd_IRA2_AH23.variant3" "Mmd_IRA3_JR11.variant3" "Mmd_IRA4_JR15.variant3" "Mmd_IRA5_JR2-F1C.variant3" "Mmd_IRA6_JR5-F1C.variant3" "Mmd_IRA7_JR7-F1C.variant3" "Mmd_IRA8_JR8-F1A.variant3" "Mmd_HEL1_HG06.variant4" "Mmd_HEL2_HG08.variant4" "Mmd_HEL3_HG13.variant4" "Mmm_CZE1_CR12.variant5" "Mmm_CZE2_CR13.variant5" "Mmm_CZE3_CR16.variant5" "Mmm_CZE4_CR23.variant5" "Mmm_CZE5_CR25.variant5" "Mmm_CZE6_CR29.variant5" "Mmm_CZE7_CR46.variant5" "Mmm_CZE8_CR59.variant5" "Mmm_KAZ1_AL1.variant6" "Mmm_KAZ2_AL16.variant6" "Mmm_KAZ3_AL19.variant6" "Mmm_KAZ4_AL33.variant6" "Mmm_KAZ5_AL38.variant6" "Mmm_KAZ6_AL40.variant6" "Mmm_KAZ7_AL41.variant6" "Mmm_KAZ8_AL42.variant6" "Mmm_AFG1_396.variant7" "Mmm_AFG2_413.variant7" "Mmm_AFG3_416.variant7" "Mmm_AFG4_424.variant7" "Mmm_AFG5_435.variant7" "Mmm_AFG6_444.variant7" "Mmc_CAST10_H36.variant8" "Mmc_CAST1_H12.variant8" "Mmc_CAST2_H14.variant8" "Mmc_CAST3_H15.variant8" "Mmc_CAST4_H24.variant8" "Mmc_CAST5_H26.variant8" "Mmc_CAST6_H27.variant8" "Mmc_CAST7_H28.variant8" "Mmc_CAST8_H30.variant8" "Mmc_CAST9_H34.variant8" "Ms_SPRE1_SP36.variant9" "Ms_SPRE2_SP39.variant9" "Ms_SPRE3_SP41.variant9" "Ms_SPRE4_SP51.variant9" "Ms_SPRE5_SP62.variant9" "Ms_SPRE6_SP68.variant9" "Ms_SPRE7_SP69.variant9" "Ms_SPRE8_SP70.variant9")

# Run bcftools consensus command on each sample
for val in ${StringArray[@]}; do
   bcftools consensus -f mm10_prdm9_reference_sequence.fa -s $val -o ${val}_prdm9_cds.fa PRDM9_CDS_AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.recode.vcf.gz
done