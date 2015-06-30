# iCAGES
This is iCAGES command line that prioritizes personalized cancer driver mutations, genes and therapies.

## Introduction
All cancers arise as a result of the acquisition of somatic mutations that drive the disease progression. It remains a challenge to identify driver mutations/genes for an individual patient and design drug therapies. To tackle this challenge, we developed iCAGES, a novel statistical framework to rapidly analyze patient-specific cancer genomic data, prioritize personalized cancer driver events and predict personalized therapies. An iCAGES web server can be downloaded from [here](http://www.github.com/WangGenomicsLab/icages-server) and installed locally.

## Download
```
wget https://github.com/WangGenomicsLab/icages/archive/v0.1.tar.gz
```

## Installation
- unzip downloaded file
```
tar -zxvf v0.1.tar.gz
```
- download and unzip database files
```
cd icages-0.1/
wget http://icages.usc.edu/download/icages/db.tar.gz
tar -zxvf db.tar.gz
```
- install necessary packages for perl
  * if you have root access, please use cpanm command to download JSON, HTTP::Request and LWP packages for perl
```
cpanm JSON
cpanm HTTP::Request
cpanm LWP
```
- install other dependencies for iCAGES 
  *  ANNOVAR
    1. please visit [ANNOVAR](http://www.openbioinformatics.org/annovar/annovar_download.html) website and download it
    2. if your current direcotry is icages-0.1, then please move annovar/ directory to ./bin diretory 
```
mv path-to-annovar/annovar/ ./bin/
```
  * DGIdb
    1. if your current directory is icages-0.1, then please create a directory under ./bin directory and name it DGIdb
    2. please visit [DGIdb](http://dgidb.genome.wustl.edu/) to read about it and download download the corresponding perl script from [here](wget https://raw.github.com/genome/dgi-db/master/files/perl_example.pl) to ./bin/DGIdb directory
```
mkdir ./bin/DGIdb
wget https://raw.github.com/genome/dgi-db/master/files/perl_example.pl -O ./bin/DGIdb/get_DrugList.pl
```
    3. please make some modifications of this get_DrugList.pl file, including these following changes:
     a. add this following line after "parse_opts();" 
``` 
open (OUT, ">$output") or die "iCAGES: cannot open file $output for writing the drugs recommended for cancer driver genes\n";
```
     b. add this following line after "'help'                  => \$help,"
```
'output:s'              => \$output 
```
     c. comment out this following line 
```
print "gene_name\tdrug_name\tinteraction_type\tsource\tgene_categories\n";
```
     d. change this following line
```
print "$gene_name\t$drug_name\t$interaction_type\t$source\t$gene_categories\n"; 
```
     into
```
print OUT "$gene_name\t$drug_name\t$interaction_type\t$source\t$gene_categories\n"; 
```
     e. change this following lines
```
print "\n" . 'Unmatched search term: ', $_->{searchTerm}, "\n";
print 'Possible suggestions: ', join(",", @{$_->{suggestions}}), "\n";
```
     into
```
print OUT "\n" . 'Unmatched search term: ', $_->{searchTerm}, "\n";
print OUT 'Possible suggestions: ', join(",", @{$_->{suggestions}}), "\n";
```
     
## Synopsis

- print help message
```
perl icages.pl --help
```

- use icages for prioritize mutations in input.txt file (in ANNOVAR input format)
```
perl icages.pl input.txt
```

- use icages for prioritize mutations in input.vcf file (in VCF input format)
```
perl icages.pl input.vcf
```

## How iCAGES works
![iCAGES pipeline](/img/icages_pipeline.png)

A flowchart to show the process of iCAGES package. The iCAGES package consists of three layers. The input file contains all variants identified from the patient; it can be either in ANNOVAR input format or in VCF format. The first layer of iCAGES prioritizes mutations. It computes three different feature scores for annotating the gene, including radial SVM score for each of its point coding mutation, CNV normalized peak score for each of its structural variation and FunSeq2 score for each of its point non-coding mutation. The second layer of iCAGES prioritizes cancer driver genes. It takes three features scores from the first layer, generates the corresponding Phenolyzer score for each mutated gene and computes an LR score for this gene (iCAGES gene score). The final level of iCAGES prioritizes targeted drugs. It first queries DGIdb for potential drugs that interact with the patient’ s mutated genes and their neighbors. Next, it calculates the joint probability for each drug being the most effective (iCAGES drug score) from three feature scores, which are iCAGES score for its direct/indirect target, normalized BioSystems probability measuring the maximum relatedness of the drugs’ direct target with each mutated gene (final target) in the patient and PubChem active probability measuring the bioactivity of the drug. The final output of iCAGES consists of three major elements, a prioritized list of mutations, a prioritized list of genes with their iCAGES gene scores as well as a prioritized list of targeted drugs with their iCAGES drug scores, all highlighted in red.

## License Agreement
By using the software, you acknowledge that you agree to the terms below:

For academic and non-profit use, you are free to fork, download, modify, distribute and use the software without restriction.

For commercial use, you are required to contact [Stevens Institute of Innovation](https://stevens.usc.edu/contact-us/) at USC directly to discuss licensing options.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Reference
Dong C, Yang H, He Z, Liu X, Wang K. **iCAGES: integrated CAncer GEnome Score for understanding personal cancer genomes**. bioRxiv doi: http://dx.doi.org/10.1101/015008

## Contact
- Coco Chengliang Dong (chenglid@usc.edu)
- Kai Wang (kaiwang@usc.edu)


