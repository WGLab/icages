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

- initialize icages (after downloaded iCAGES, please first INITIALIZE it before using it): 
```
perl icagesInitiate.pl
```

- use icages for prioritize mutations in input.txt file (in ANNOVAR input format)
```
perl icages.pl input.txt
```

- use icages for prioritize mutations in input.vcf file (in VCF input format)
```
perl icages.pl input.vcf
```


## License Agreement
By using the software, you acknowledge that you agree to the terms below:

For academic and non-profit use, you are free to fork, download, modify, distribute and use the software without restriction.

For commercial use, you are required to contact [Stevens Institute of Innovation](https://stevens.usc.edu/contact-us/) at USC directly to discuss licensing options.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Contact
- Coco Chengliang Dong (chenglid@usc.edu)
- Kai Wang (kaiwang@usc.edu)


