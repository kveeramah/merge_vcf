# merge_vcf
This program merges individual emit all vcf files created by the aDNA_GenoCaller and GenoCaller_indent scripts for biallelic SNPs

pysam needs to be installed. Two input files are needed.

The first input file needs to be a tab delimited bed file containing the list of SNPs to be included in the output vcf, one per line, of the form <chrom><start><end><all1><all2><snp_name>

e.g. 2\t136608645\t136608646\tG\tA\t2:136608646\n

One of two alleles should be the reference allele, though the order does not matter

The second input file contains a list of individual emit_all vcf files to be merged (exact location with respect to working directory), 

followed by the assigned sample name, tab seperated.

e.g.
../2072_LCT.RG.LCT.bed.aDNA.emit_all.vcf\tsamp_2072
../84001_LCT.RG.LCT.bed.aDNA.emit_all.vcf\tsamp_84001
../84005_LCT.RG.LCT.bed.aDNA.emit_all.vcf\tsamp_84005

to run type:
merge_vcf.py <SNP_bedfile> list_of_emit_alls> <reference genome>

Output will be a merged vcf.gz file.
