## ANNOVAR input with new prefix and new direcotry for output

To change prefix, please use option `-p` or `--prefix` and to change directory where your output will be generated, please use option `--outputdir`. Note that if you have other forms of input, such as VCF format and BED format, the syntax is the same.

```
[cocodong@biocluster ~/]$ head input.txt 
1	12919840	12919840	T	C
1	35332717	35332717	C	A
1	55148456	55148456	G	T
1	70504789	70504789	C	T
1	167059520	167059520	A	T
1	182496864	182496864	A	T
1	197073351	197073351	C	T
1	216373211	216373211	G	T
10	37490170	37490170	G	A
10	56089432	56089432	A	C
[cocodong@biocluster ~/]$ icages.pl input.txt -p newname --outputdir newoutputdir

```

## ANNOVAR input annotated with hg38

To change database version, please use option `--buildver`. Note that if you have other forms of input, such as VCF format and BED format, the syntax is the same.

```
[cocodong@biocluster ~/]$ head input.txt 
1	12919840	12919840	T	C
1	35332717	35332717	C	A
1	55148456	55148456	G	T
1	70504789	70504789	C	T
1	167059520	167059520	A	T
1	182496864	182496864	A	T
1	197073351	197073351	C	T
1	216373211	216373211	G	T
10	37490170	37490170	G	A
10	56089432	56089432	A	C
[cocodong@biocluster ~/]$ icages.pl input.txt --buildver hg38

```

## VCF input with one sample which contains both germline mutations and mutations in his/her tumor

If you do not have somatic mutations for one sample in VCF file, but what you have is a VCF file that contains both germline mutations and mutations in cancer for this sample, then you can specify the headers for germline mutations using options `-t` or `--tumor` and specify the headers for tumor mutations using options `-g` or `--germline`. iCAGES will be able to extract somatic mutations from this VCF file and carry on downstream analysis for you. In this example, the input file is a VCF file that contains tumor mutations with header "tumor" and germline mutations with header "germline", all annotated with reference genome version of hg19.
```
[cocodong@biocluster ~/]$ cat input.vcf
##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tumor	germline
1	12919840	.	T	C	.	.	.	GT	1|1	0|0
1	35332717	.	C	A	.	.	.	GT	1|1	0|0
1	55148456	.	G	T	.	.	.	GT	1|1	0|0
1	70504789	.	C	T	.	.	.	GT	1|1	0|0
1	167059520	.	A	T	.	.	.	GT	1|1	0|0
1	182496864	.	A	T	.	.	.	GT	1|1	0|0
1	197073351	.	C	T	.	.	.	GT	1|1	0|0
1	216373211	.	G	T	.	.	.	GT	1|1	0|0
10	37490170	.	G	A	.	.	.	GT	1|1	0|0
10	56089432	.	A	C	.	.	.	GT	1|1	0|0
...
[cocodong@biocluster ~/]$ icages.pl input.vcf -t tumor -g germline
```

## VCF input with multiple samples which contains both germline mutations and tumor mutations

iCAGES is a personalized cancer driver analysis pipeline, so it only does analysis for ONE single patient. But if what you have is a VCF file that contains both germline mutations and tumor mutations for multiple individuals, then you can specify the headers for germline mutations for the patient of your interest using options `-t` or `--tumor` and specify the headers for tumor mutations for the patient of your interest using options `-g` or `--germline`. iCAGES will be able to extract somatic mutations for this particular individual from this VCF file and carry on downstream analysis for you. In this example, the input file is a VCF file that contains mutations from two individuals Sapmle1 and Sample2, each of them have both tumor mutations and germline mutations with slightly different headers, all annotated with reference genome version of hg19. By specifying headers for Sample1, iCAGES analyzes this sample for you.
```
[cocodong@biocluster ~/]$ cat input.vcf
##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1Tumor	Sample1Germline	Sample2Tumor	Sample2Germline
1	12919840	.	T	C	.	.	.	GT	1|1	0|0	1|1	0|0
1	35332717	.	C	A	.	.	.	GT	1|1	0|0	1|1	0|0
1	55148456	.	G	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	70504789	.	C	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	167059520	.	A	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	182496864	.	A	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	197073351	.	C	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	216373211	.	G	T	.	.	.	GT	1|1	0|0	1|1	0|0
10	37490170	.	G	A	.	.	.	GT	1|1	0|0	1|1	0|0
10	56089432	.	A	C	.	.	.	GT	1|1	0|0	1|1	0|0
...
[cocodong@biocluster ~/]$ icages.pl input.vcf -t Sample1Tumor -g Sample1Germline
```


## VCF input with multiple samples which contains both germline mutations and tumor mutations

iCAGES is a personalized cancer driver analysis pipeline, so it only does analysis for ONE single patient. But if what you have is a VCF file that contains both germline mutations and tumor mutations for multiple individuals, then you can specify the headers for germline mutations for the patient of your interest using options `-t` or `--tumor` and specify the headers for tumor mutations for the patient of your interest using options `-g` or `--germline`. iCAGES will be able to extract somatic mutations for this particular individual from this VCF file and carry on downstream analysis for you. In this example, the input file is a VCF file that contains mutations from two individuals Sapmle1 and Sample2, each of them have both tumor mutations and germline mutations with slightly different headers, all annotated with reference genome version of hg19. By specifying headers for Sample1, iCAGES analyzes this sample for you.
```
[cocodong@biocluster ~/]$ cat input.vcf
##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1Tumor	Sample1Germline	Sample2Tumor	Sample2Germline
1	12919840	.	T	C	.	.	.	GT	1|1	0|0	1|1	0|0
1	35332717	.	C	A	.	.	.	GT	1|1	0|0	1|1	0|0
1	55148456	.	G	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	70504789	.	C	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	167059520	.	A	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	182496864	.	A	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	197073351	.	C	T	.	.	.	GT	1|1	0|0	1|1	0|0
1	216373211	.	G	T	.	.	.	GT	1|1	0|0	1|1	0|0
10	37490170	.	G	A	.	.	.	GT	1|1	0|0	1|1	0|0
10	56089432	.	A	C	.	.	.	GT	1|1	0|0	1|1	0|0
...
[cocodong@biocluster ~/]$ icages.pl input.vcf -t Sample1Tumor -g Sample1Germline
```

## VCF input with multiple samples which contains only somatic mutations

Again, iCAGES is a personalized cancer driver analysis pipeline, so it only does analysis for ONE single patient. But if what you have is a VCF file that contains somatic mutations for multiple individuals, then you can specify the header for the patient of your interest using options `-i` or `--id`. iCAGES will be able to extract somatic mutations for this particular individual from this VCF file and carry on downstream analysis for you. In this example, the input file is a VCF file that contains somatic mutations from two individuals Sapmle1 and Sample2, all annotated with reference genome version of hg19. By specifying header for Sample1, iCAGES analyzes this sample for you.
```
[cocodong@biocluster ~/]$ cat input.vcf
##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
1	12919840	.	T	C	.	.	.	GT	1|1	0|0
1	35332717	.	C	A	.	.	.	GT	1|1	1|1
1	55148456	.	G	T	.	.	.	GT	1|1	0|0
1	70504789	.	C	T	.	.	.	GT	1|1	1|1
1	167059520	.	A	T	.	.	.	GT	1|1	0|0
1	182496864	.	A	T	.	.	.	GT	0|0	0|0
1	197073351	.	C	T	.	.	.	GT	1|1	1|1
1	216373211	.	G	T	.	.	.	GT	1|1	0|0
10	37490170	.	G	A	.	.	.	GT	1|1	0|0
10	56089432	.	A	C	.	.	.	GT	1|1	0|0
...
[cocodong@biocluster ~/]$ icages.pl input.vcf -i Sample1
```

## VCF input with multiple samples and BED files with additional structural variations

VCF has immaure development of annotation on structural variations. In order to better annotate personal cancer mutation profiles, we made iCAGES to support additional BED file input, which profiles structural variations, using options `-b` or `--bed`. iCAGES will be able to combine information from VCF files and BED files to do downstream data analysis for you. In this example, the input files are a VCF file that contains somatic mutations from two individuals Sapmle1 and Sample2, all annotated with reference genome version of hg19 and a BED file that contains coordinates of structural varations. This exemplary BED file is also provided in the package.
```
[cocodong@biocluster ~/]$ cat input.vcf
##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
1	12919840	.	T	C	.	.	.	GT	1|1	0|0
1	35332717	.	C	A	.	.	.	GT	1|1	1|1
1	55148456	.	G	T	.	.	.	GT	1|1	0|0
1	70504789	.	C	T	.	.	.	GT	1|1	1|1
1	167059520	.	A	T	.	.	.	GT	1|1	0|0
1	182496864	.	A	T	.	.	.	GT	0|0	0|0
1	197073351	.	C	T	.	.	.	GT	1|1	1|1
1	216373211	.	G	T	.	.	.	GT	1|1	0|0
10	37490170	.	G	A	.	.	.	GT	1|1	0|0
10	56089432	.	A	C	.	.	.	GT	1|1	0|0
...
[cocodong@biocluster ~/]$ cat input.bed
chr12	85865797	85887628
chr20	15052592	15071191
chr16	87340388	87349798
chr2		213000509	213007522
[cocodong@biocluster ~/]$ icages.pl input.vcf -i Sample1 -b input.bed
```




---

<div id="disqus_thread"></div>
<script type="text/javascript">
/* * * CONFIGURATION VARIABLES * * */
var disqus_shortname = 'icages';
var disqus_identifier = 'example';
var disquss_title = 'iCAGES Examples';

/* * * DON'T EDIT BELOW THIS LINE * * */
(function() {
var dsq = document.createElement('script'); dsq.type = 'text/javascript'; dsq.async = true;
dsq.src = '//' + disqus_shortname + '.disqus.com/embed.js';
(document.getElementsByTagName('head')[0] || document.getElementsByTagName('body')[0]).appendChild(dsq);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript" rel="nofollow">comments powered by Disqus.</a></noscript>


