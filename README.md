# iCAGES
This is iCAGES command line that prioritizes personalized cancer driver mutations, genes and therapies.

## Introduction
All cancers arise as a result of the acquisition of somatic mutations that drive the disease progression. It remains a challenge to identify driver mutations/genes for an individual patient and design drug therapies. To tackle this challenge, we developed iCAGES, a novel statistical framework to rapidly analyze patient-specific cancer genomic data, prioritize personalized cancer driver events and predict personalized therapies.

## Synopsis

- download iCAGES
```
git clone https://github.com/WangGenomicsLab/icages.git
```

- print help message
```
perl icages.pl --help
```

- initialize icages (after downloaded iCAGES, please first INITIALIZE it before using it). Please read the documentation on downloading databases.

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

## Contact
- Coco Chengliang Dong (chenglid@usc.edu)
- Kai Wang (kaiwang@usc.edu)


