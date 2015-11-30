## iCAGES main package

Please join the iCAGES mailing list at google groups [here](https://groups.google.com/forum/?hl=en#!forum/icages) to receive announcements on software updates.

The latest version of iCAGES (2015Nov27) can be downloaded [here](https://github.com/WangGenomicsLab/icages/releases/tag/v1.0.1).

iCAGES is written in Perl and can be run as a standalone application on diverse hardware systems where standard Perl modules are installed.

## Download

```
wget https://codeload.github.com/WangGenomicsLab/icages/tar.gz/(version)
```

## Installation

- Unzip downloaded file

```
tar -zxvf icages-(version).tar.gz
mv icages-(version) icages
```

- Download and unzip database files

```
cd icages/
wget http://icages.usc.edu/download/icages/db.tar.gz
tar -zxvf db.tar.gz
```

- Install necessary packages for perl. If you have root access, please use cpanm command to download JSON, HTTP::Request and LWP packages for perl

```
cpanm JSON
cpanm HTTP::Request
cpanm LWP
```

- Install the first dependency for iCAGES, ANNOVAR. Please visit [ANNOVAR](http://www.openbioinformatics.org/annovar/annovar_download.html) website and download it. If your current direcotry is icages, then please move annovar/ directory to ./bin diretory 
```
mv path-to-annovar/annovar/ ./bin/
```

- Install the second dependency for iCAGES, DGIdb. If your current directory is icages, then please create a directory under ./bin directory and name it DGIdb.Please visit [DGIdb](http://dgidb.genome.wustl.edu/) to read about it and download download the corresponding perl script from [here](wget https://raw.github.com/genome/dgi-db/master/files/perl_example.pl) to ./bin/DGIdb directory

```
mkdir ./bin/DGIdb
wget https://raw.github.com/genome/dgi-db/master/files/perl_example.pl -O ./bin/DGIdb/get_DrugList.pl
```

- Please make some modifications of this get_DrugList.pl file. First, add this following line after "parse_opts();" 

``` 
open (OUT, ">$output") or die "iCAGES: cannot open file $output for writing the drugs recommended for cancer driver genes\n";
```

- Then, add this following line after "'help' => \$help,"

```
'output:s'    => \$output 
```

- Next, comment out this following line 

```
print "gene_name\tdrug_name\tinteraction_type\tsource\tgene_categories\n";
```
into
```
# print "gene_name\tdrug_name\tinteraction_type\tsource\tgene_categories\n";
```

- Next, change this following line

```
print "$gene_name\t$drug_name\t$interaction_type\t$source\t$gene_categories\n"; 
```
into
```
print OUT "$gene_name\t$drug_name\t$interaction_type\t$source\t$gene_categories\n"; 
```

- And change this following lines

```
print "\n" . 'Unmatched search term: ', $_->{searchTerm}, "\n";
print 'Possible suggestions: ', join(",", @{$_->{suggestions}}), "\n";
```
into 
```
print OUT "\n" . 'Unmatched search term: ', $_->{searchTerm}, "\n";
print OUT 'Possible suggestions: ', join(",", @{$_->{suggestions}}), "\n";
```

- Install the third dependency for iCAGES, vcftools. Asuming you are already in icages-(version)/bin/ directory, download vcftools through sourceforge

```
wget http://iweb.dl.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz
tar -zxvf vcftools_0.1.12b.tar.gz 
mv vcftools_0.1.12b/ vcftools/
rm vcftools_0.1.12b.tar.gz
```

- Install vcftools by typing the following command

```
cd vcftools
make
```

- Install the fourth dependency for iCAGES, bedtools. Asuming you are already in icages-(version)/bin/ directory,

```
wget https://codeload.github.com/arq5x/bedtools2/tar.gz/v2.25.0
tar -zxvf v2.2.25.0.tar.gz
mv v2.2.25.0 bedtools
rm v2.2.25.0.tar.gz
cd bedtools
make
```

## Additional databases

Initial databases for iCAGES only includes hg19 reference genome for human. In order to annotate variants with hg18 or hg38 reference genomes, please download these additional databases compiled for these two versions of references.

- hg18

```
cd icages/db/
wget http://icages.usc.edu/download/icages/db_hg18.tar.gz
tar -zxvf db_hg18.tar.gz
```

- hg38

```
cd icages/db/
wget http://icages.usc.edu/download/icages/db_hg38.tar.gz
tar -zxvf db_hg18.tar.gz
```



---

<div id="disqus_thread"></div>
<script type="text/javascript">
/* * * CONFIGURATION VARIABLES * * */
var disqus_shortname = 'icages';
var disqus_identifier = 'download';
var disquss_title = 'iCAGES Download';

/* * * DON'T EDIT BELOW THIS LINE * * */
(function() {
var dsq = document.createElement('script'); dsq.type = 'text/javascript'; dsq.async = true;
dsq.src = '//' + disqus_shortname + '.disqus.com/embed.js';
(document.getElementsByTagName('head')[0] || document.getElementsByTagName('body')[0]).appendChild(dsq);
})();
</script>
<noscript>Please enable JavaScript to view the <a href="https://disqus.com/?ref_noscript" rel="nofollow">comments powered by Disqus.</a></noscript>

