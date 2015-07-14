## ANNOVAR input 

For beginners, the easiest way to use iCAGES is to annotate somatic mutations in [ANNOVAR](http://annovar.openbioinformatics.org/) input format with reference genome version hg19. This exemplary input file is provided in iCAGES package.
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
[cocodong@biocluster ~/]$ icages.pl input.txt

```

## VCF input with one sample with his/her somatic mutations only

If you have somatic mutations for one sample in VCF file format with reference genome version hg19, iCAGES can automaticaly detect the input format and analyze your data. This exemplary input file is also provided in iCAGES package.
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
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  tumor
1       12919840        .       T       C       .       .       .       GT      1|1
1       35332717        .       C       A       .       .       .       GT      1|1
1       55148456        .       G       T       .       .       .       GT      1|1
1       70504789        .       C       T       .       .       .       GT      1|1
1       167059520       .       A       T       .       .       .       GT      1|1
1       182496864       .       A       T       .       .       .       GT      1|1
1       197073351       .       C       T       .       .       .       GT      1|1
1       216373211       .       G       T       .       .       .       GT      1|1
10      37490170        .       G       A       .       .       .       GT      1|1
10      56089432        .       A       C       .       .       .       GT      1|1
...
[cocodong@biocluster ~/]$ icages.pl input.vcf
```

## BED input with one sample with his/her somatic structural variations only

If you only have BED files for all structural variations detected in this patient with reference genome version hg19, iCAGES can also automatically detect the input format of your data and carry on downstream analysis. This exemplary input file is also provided in iCAGES package.

```
[cocodong@biocluster ~/]$ head input.bed
chr12	85865797	85887628
chr20	15052592	15071191
chr16	87340388	87349798
chr2		213000509	213007522
[cocodong@biocluster ~/]$ icages.pl input.bed
```


---

<div id="disqus_thread"></div>
<script type="text/javascript">
/* * * CONFIGURATION VARIABLES * * */
var disqus_shortname = 'icages';

/* * * DON'T EDIT BELOW THIS LINE * * */
(function() {
var dsq = document.createElement('script'); dsq.type = 'text/javascript'; dsq.async = true;
dsq.src = '//' + disqus_shortname + '.disqus.com/embed.js';
(document.getElementsByTagName('head')[0] || document.getElementsByTagName('body')[0]).appendChild(dsq);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript" rel="nofollow">comments powered by Disqus.</a></noscript>
