#!/bin/bash

## 12th November 2020
wget http://berglandlab.uvadcos.io/vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz

# 1st December 2020
wget http://berglandlab.uvadcos.io/vcf/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz


## Remove problematic populations
grep -v "##" dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz > nohead_dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
awk awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$12"\t"$13"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36"\t"$38"\t"$39"\t"$40"\t"$41"\t"$42"\t"$43"\t"$44"\t"$45"\t"$46"\t"$47"\t"$48"\t"$49"\t"$50"\t"$51"\t"$52"\t"$53"\t"$54"\t"$55"\t"$56"\t"$57"\t"$58"\t"$59"\t"$60"\t"$61"\t"$62"\t"$63"\t"$64"\t"$65"\t"$66"\t"$67"\t"$68"\t"$69"\t"$70"\t"$71"\t"$72"\t"$73"\t"$74"\t"$75"\t"$76"\t"$77"\t"$78"\t"$79"\t"$80"\t"$81"\t"$83"\t"$85"\t"$86"\t"$87"\t"$88"\t"$91"\t"$92"\t"$93"\t"$94"\t"$95"\t"$96"\t"$97"\t"$98"\t"$99"\t"$100"\t"$101"\t"$102"\t"$103"\t"$104"\t"$105"\t"$106"\t"$107"\t"$108"\t"$109"\t"$110"\t"$112"\t"$113"\t"$114"\t"$115"\t"$116"\t"$117"\t"$118"\t"$119"\t"$120"\t"$121"\t"$122"\t"$124"\t"$125"\t"$126"\t"$127"\t"$128"\t"$129"\t"$130"\t"$131"\t"$132"\t"$133"\t"$134"\t"$135"\t"$136"\t"$137"\t"$138"\t"$139"\t"$140"\t"$141"\t"$142"\t"$143"\t"$144"\t"$145"\t"$146"\t"$147"\t"$148"\t"$150"\t"$151"\t"$152"\t"$153"\t"$154"\t"$155"\t"$156"\t"$157"\t"$158"\t"$159"\t"$162"\t"$163"\t"$164"\t"$165"\t"$166"\t"$167"\t"$168"\t"$169"\t"$170"\t"$171"\t"$172"\t"$173"\t"$174"\t"$175"\t"$178"\t"$179"\t"$180"\t"$181"\t"$182"\t"$183"\t"$184"\t"$185"\t"$186"\t"$187"\t"$188"\t"$189"\t"$190"\t"$191"\t"$192"\t"$193"\t"$196"\t"$197"\t"$198"\t"$199"\t"$201"\t"$202"\t"$203"\t"$204"\t"$205"\t"$206"\t"$207"\t"$208"\t"$209"\t"$210"\t"$211"\t"$212"\t"$213"\t"$214"\t"$215"\t"$216"\t"$217"\t"$218"\t"$220"\t"$221"\t"$222"\t"$223"\t"$224"\t"$225"\t"$226"\t"$227"\t"$228"\t"$229"\t"$230"\t"$231"\t"$232"\t"$233"\t"$234"\t"$235"\t"$236"\t"$237"\t"$238"\t"$239"\t"$240"\t"$242"\t"$243"\t"$244"\t"$245"\t"$246"\t"$247"\t"$248"\t"$249"\t"$250"\t"$251"\t"$252"\t"$253"\t"$254"\t"$255}' nohead_dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf | gzip > nohead_dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz


## MAF filter
module load Python/3.8.2-GCCcore-9.3.0

# Obtain all positions including monomorphic for all populations
python3  FilterSNAPEByMAF_MB.py \
--VCF  dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz \
--MAF 0.05 \
| gzip > MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz

# Obtain the filtered vcf without the problematic populations to have the final SNP dataset needed
gunzip nohead_dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz
cat header_dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf nohead_dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf | gzip > dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz
gzip nohead_dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf

python3 FilterSNAPEByMAF.py \
--VCF  dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz \
--MAF 0.05 \
| gzip > MAF005.dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz

## Get positions
zcat MAF005.dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz|  awk '{print $1$2}' > POSITIONS_MAF005.dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
zcat MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz | awk '{print $0"\t"$1$2}' > withPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf

awk -F '\t' 'NR==FNR {id[$1]; next} $256 in id' POSITIONS_MAF005.dest_nonproblematic.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf  withPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf | gzip > filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz


## Splitting SNAPE file into chromosomes
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz| awk '$1=="2L" {print $0}' > 2L_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz| awk '$1=="2R" {print $0}' > 2R_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz| awk '$1=="3L" {print $0}' > 3L_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz| awk '$1=="3R" {print $0}' > 3R_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz| awk '$1=="4" {print $0}' > 4_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz| awk '$1=="Y" {print $0}' > Y_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz| awk '$1=="X" {print $0}' > X_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf

## Splitting PoolSNP file into chromosomes
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | awk '$1=="2L" {print $0}' > 2L_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | awk '$1=="2R" {print $0}' > 2R_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | awk '$1=="3L" {print $0}' > 3L_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | awk '$1=="3R" {print $0}' > 3R_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | awk '$1=="4" {print $0}' > 4_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | awk '$1=="Y" {print $0}' > Y_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | awk '$1=="X" {print $0}' > X_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf

## Headers
zcat  dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | grep "#CHROM" > header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf
zcat  filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz | grep "#CHROM" > header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf


module load Python/2.7.15-foss-2018b
## 2L parsing
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 2L_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 2L_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf > 2L_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 2L_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 2L_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf > 2L_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf
cat header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 2L_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf > 2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf
cat header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 2L_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf > 2L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf
paste 2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > 2L_PoolSNP_SNAPE_001.50.ann.vcf
# Per population
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf > ${i}_2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf; done
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 2L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_2L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf; done
# Create the 2-2 data sets
for i in {10..255}; do paste 1_2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2_2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_2L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf ${i}_2L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf ${i}_2L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_2L_header_PoolSNP_SNAPE.vcf ; done
for i in {10..255}; do python caller_comparison_per_pop.py ${i}_2L_header_PoolSNP_SNAPE.vcf ${i}_2L_header_PoolSNP_SNAPE.stats.txt; done
for i in {10..255}; do tail -n1 ${i}_2L_header_PoolSNP_SNAPE.stats.txt | awk -v V=$i '{print "'${i}'\t"$0}'; done > all_pops_2L_header_PoolSNP_SNAPE.stats.txt


## 2R parsing
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 2R_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 2R_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf > 2R_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 2R_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 2R_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf > 2R_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf
cat header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 2R_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf > 2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf
cat header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 2R_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf > 2R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf
paste 2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > 2R_PoolSNP_SNAPE_001.50.ann.vcf
# Per population
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf > ${i}_2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf; done
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 2R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_2R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf; done
# Create the 2-2 data sets
for i in {10..255}; do paste 1_2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2_2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_2R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf ${i}_2R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf ${i}_2R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_2R_header_PoolSNP_SNAPE.vcf; done
for i in {10..255}; do python caller_comparison_per_pop.py ${i}_2R_header_PoolSNP_SNAPE.vcf ${i}_2R_header_PoolSNP_SNAPE.stats.txt; done
for i in {10..255}; do tail -n1 ${i}_2R_header_PoolSNP_SNAPE.stats.txt | awk -v V=$i '{print "'${i}'\t"$0}'; done > all_pops_2R_header_PoolSNP_SNAPE.stats.txt


## 3L parsing
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 3L_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 3L_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf > 3L_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 3L_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 3L_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf > 3L_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf
cat header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 3L_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf > 3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf
cat header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 3L_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf > 3L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf
paste 3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 3L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > 3L_PoolSNP_SNAPE_001.50.ann.vcf
# Per population
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf > ${i}_3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf; done
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 3L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_3L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf; done
# Create the 2-2 data sets
for i in {10..255}; do paste 1_3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2_3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_3L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf ${i}_3L_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf ${i}_3L_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_3L_header_PoolSNP_SNAPE.vcf ; done
for i in {10..255}; do python caller_comparison_per_pop.py ${i}_3L_header_PoolSNP_SNAPE.vcf ${i}_3L_header_PoolSNP_SNAPE.stats.txt; done
for i in {10..255}; do tail -n1 ${i}_3L_header_PoolSNP_SNAPE.stats.txt | awk -v V=$i '{print "'${i}'\t"$0}'; done > all_pops_3L_header_PoolSNP_SNAPE.stats.txt


## 3R parsing
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 3R_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 3R_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf > 3R_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 3R_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 3R_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf > 3R_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf
cat header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 3R_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf > 3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf
cat header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 3R_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf > 3R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf
paste 3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 3R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > 3R_PoolSNP_SNAPE_001.50.ann.vcf
# Per population
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf > ${i}_3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf; done
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 3R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_3R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf; done
# Create the 2-2 data sets
for i in {10..255}; do paste 1_3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2_3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_3R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf ${i}_3R_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf ${i}_3R_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_3R_header_PoolSNP_SNAPE.vcf ; done
for i in {10..255}; do python caller_comparison_per_pop.py ${i}_3R_header_PoolSNP_SNAPE.vcf ${i}_3R_header_PoolSNP_SNAPE.stats.txt; done
for i in {10..255}; do tail -n1 ${i}_3R_header_PoolSNP_SNAPE.stats.txt | awk -v V=$i '{print "'${i}'\t"$0}'; done > all_pops_3R_header_PoolSNP_SNAPE.stats.txt


## 4 parsing
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 4_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 4_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf > 4_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' 4_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 4_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf > 4_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf
cat header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf 4_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf > 4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf
cat header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf 4_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf > 4_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf
paste 4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > 4_PoolSNP_SNAPE_001.50.ann.vcf
# Per population
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf > ${i}_4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf; done
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' 4_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_4_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf; done
# Create the 2-2 data sets
for i in {10..255}; do paste 1_4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2_4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_4_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf ${i}_4_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf ${i}_4_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_4_header_PoolSNP_SNAPE.vcf ; done
for i in {10..255}; do python caller_comparison_per_pop.py ${i}_4_header_PoolSNP_SNAPE.vcf ${i}_4_header_PoolSNP_SNAPE.stats.txt; done
for i in {10..255}; do tail -n1 ${i}_4_header_PoolSNP_SNAPE.stats.txt | awk -v V=$i '{print "'${i}'\t"$0}'; done > all_pops_4_header_PoolSNP_SNAPE.stats.txt


## Y parsing
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' Y_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf Y_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf > Y_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' Y_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf Y_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf > Y_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf
cat header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf Y_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf > Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf
cat header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf Y_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf > Y_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf
paste Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf Y_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > Y_PoolSNP_SNAPE_001.50.ann.vcf
# Per population
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf > ${i}_Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf; done
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' Y_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_Y_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf; done
# Create the 2-2 data sets
for i in {10..255}; do paste 1_Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2_Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_Y_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf ${i}_Y_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf ${i}_Y_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_Y_header_PoolSNP_SNAPE.vcf ; done
for i in {10..255}; do python caller_comparison_per_pop.py ${i}_Y_header_PoolSNP_SNAPE.vcf ${i}_Y_header_PoolSNP_SNAPE.stats.txt; done
for i in {10..255}; do tail -n1 ${i}_Y_header_PoolSNP_SNAPE.stats.txt | awk -v V=$i '{print "'${i}'\t"$0}'; done > all_pops_Y_header_PoolSNP_SNAPE.stats.txt


## X parsing
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' X_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf X_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf > X_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf
awk -F '\t' 'NR==FNR {id[$2]; next} $2 in id' X_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf X_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf > X_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf
cat header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf X_positions_dest.PoolSeq.filteredPoolSNP.001.50.ann.vcf > X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf
cat header_filteredPOS.MAF005.MB.dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf X_positions_dest.PoolSeq.filteredSNAPE.001.50.ann.vcf > X_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf
paste X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf X_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > X_PoolSNP_SNAPE_001.50.ann.vcf
# Per population
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf > ${i}_X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf; done
for i in $(seq 255 $END); do awk -v V=$i '{print $V}' X_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_X_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf; done
# Create the 2-2 data sets
for i in {10..255}; do paste 1_X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 2_X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 4_X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf 5_X_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf ${i}_X_header_positions_dest.PoolSeq.PoolSNP.001.50.ann.vcf ${i}_X_header_positions_dest.PoolSeq.SNAPE.001.50.ann.vcf > ${i}_X_header_PoolSNP_SNAPE.vcf ; done
for i in {10..255}; do python caller_comparison_per_pop.py ${i}_X_header_PoolSNP_SNAPE.vcf ${i}_X_header_PoolSNP_SNAPE.stats.txt; done
for i in {10..255}; do tail -n1 ${i}_X_header_PoolSNP_SNAPE.stats.txt | awk -v V=$i '{print "'${i}'\t"$0}'; done > all_pops_X_header_PoolSNP_SNAPE.stats.txt

