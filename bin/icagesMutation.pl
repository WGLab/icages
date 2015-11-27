#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);	
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my ($annovarInputFile,$rawInputFile, $inputDir ,$icagesLocation , $tumor ,$germline, $id ,$prefix, $bed, $hg, $expression);
my $nowString;
my (%sup, %onc);

######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################
$rawInputFile = $ARGV[0];
$inputDir = $ARGV[1];
$icagesLocation = $ARGV[2];
$tumor = $ARGV[3];
$germline = $ARGV[4];
$id = $ARGV[5];
$prefix = $ARGV[6];
$bed = $ARGV[7];
$hg = $ARGV[8];
$expression = $ARGV[9];
$nowString = localtime();
$annovarInputFile = &runAnnovar($rawInputFile, $inputDir ,$icagesLocation ,$tumor ,$germline ,$id, $prefix, $bed, $hg , $expression);
&processAnnovar($annovarInputFile, $hg, $icagesLocation);

######################################################################################################################################
############################################################# subroutines ############################################################
######################################################################################################################################

sub runAnnovar {
    print "NOTICE: start runing iCAGES packge at $nowString\n";
    my ($rawInputFile, $annovarInputFile);                                             #ANNOVAR input files
    my ($icagesLocation, $callAnnotateVariation, $DBLocation);                         #ANNOVAR commands
    my ($tumor, $germline ,$id ,$prefix);                                     #VCF conversion paramters & prefix for output
    my ($radialSVMDB, $radialSVMIndex, $funseq2DB, $funseq2Index, $refGeneDB, $refGeneIndex, $cnvDB);           #ANNOVAR DB files: iCAGES score (index), refGene (fasta), dbSNP
    my ($log, $annovarLog);
    my $nowString = localtime;                                                                                  #SYSTEM local time
    my $bed;    # location of bed file
    my $expression; # location of expression log change file
    $rawInputFile = shift;
    $inputDir = shift;
    $icagesLocation = shift;
    $tumor = shift;
    $germline = shift;
    $id = shift;
    $prefix = shift;
    $bed = shift;
    $hg = shift;
    $expression = shift;
    $callAnnotateVariation = $icagesLocation . "bin/annovar/annotate_variation.pl";
    $annovarInputFile = $inputDir . "/" . $prefix . ".annovar";
    $DBLocation = $icagesLocation . "db/";
    &formatConvert($rawInputFile, $annovarInputFile, $icagesLocation, $tumor, $germline , $id , $bed , $expression);
    &divideMutation($annovarInputFile);
    &loadDatabase($DBLocation);
    &annotateMutation($icagesLocation, $annovarInputFile, $hg);
    return $annovarInputFile;
}


sub processAnnovar{
    print "NOTICE: start processing output from ANNOVAR\n";
    my $annovarInputFile = shift;
    my $hg = shift;
    my $icagesLocation = shift;
    # still have to load onc gene set and suppressor gene set
    my $oncDB = $icagesLocation . "/db/oncogene.gene";
    my $supDB = $icagesLocation . "/db/suppressor.gene";
    my (%onc, %sup);
    open(ONCDB, "$oncDB") or die "ERROR: cannot open $oncDB\n";
    open(SUPDB, "$supDB") or die "ERROR: cannot open $supDB\n";
    while(<ONCDB>){
	chomp;
	$onc{$_} =1;
    }
    close ONCDB;
    while(<SUPDB>){
	chomp;
	$sup{$_} = 1;
    }
    close SUPDB;
    my $annovarVariantFunction = $annovarInputFile . ".variant_function";
    # create a temp file to store bed file of variant function : chr start end score
    my $genebed = $annovarInputFile . ".variant_function.bed";
    open(GENEFORBED, "$annovarVariantFunction") or die "ERROR: cannot open $annovarVariantFunction\n";
    open(GENEBED, ">$genebed") or die "ERROR: cannot create $genebed file\n";
    my $lastline = "" ;
    while(<GENEFORBED>){
	chomp;
	my @line = split("\t", $_);
	my $printout = "$line[2]\t$line[3]\t$line[4]"; 
	if($printout eq $lastline){
	    next;
	}else{
	    print GENEBED "$line[2]\t$line[3]\t$line[4]\n";
	}
	$lastline = $printout;
    }
    close GENEBED;
    close GENEFORBED;
    my $annovarExonVariantFunction = $annovarInputFile . ".exonic_variant_function";
    my $annovarRadialSVM = $annovarInputFile . ".snp." . $hg . "_iCAGES_dropped";
    my $annovarCNV = $annovarInputFile . ".cnv." . $hg . "_cnv";
    # create a temp file to store bed file of cnv with this format : chr start end score
    my $cnvbed = $annovarInputFile . ".cnv." . $hg . "_cnv.bed";

    # create a final file to store the final result of bedtools intersect
    my $cnvfinal = $annovarInputFile . ".cnv.final";

    my $annovarFunseq2 = $annovarInputFile . ".snp." . $hg . "_funseq2_dropped";
    my $icagesMutations = $annovarInputFile . ".icagesMutations.csv";
    # add bedtools 
    my $bedtools = $icagesLocation . "/bin/bedtools/bin/bedtools";
    open(GENE, "$annovarVariantFunction") or die "ERROR: cannot open file $annovarVariantFunction\n";
    open(EXON, "$annovarExonVariantFunction") or die "ERROR: cannot open file $annovarExonVariantFunction\n";
    if(!-e $annovarCNV){
	!system("touch $annovarCNV") or die "ERROR: cannot create file $annovarCNV\n";
    }
    if(!-e $annovarRadialSVM){
	!system("touch $annovarRadialSVM") or die "ERROR: cannot create file $annovarRadialSVM\n";
    }
    if(!-e $annovarFunseq2){
	!system("touch $annovarFunseq2") or die "ERROR: cannot create file $annovarFunseq2\n";
    }

    open(CNV, "$annovarCNV") or die "ERROR: cannot open file $annovarCNV\n";
    open(RADIAL, "$annovarRadialSVM") or die "ERROR: cannot open file $annovarRadialSVM\n";
    open(FUNSEQ, "$annovarFunseq2") or die "ERROR: cannot open file $annovarFunseq2\n";
    open(OUT, ">$icagesMutations") or die "ERROR: cannot open file $icagesMutations\n";
    my (%radialSVM, %funseq, %cnv, %exon);
    my (%pointcoding);
    my %icagesMutations;
    
    ######## count location information
    my $exonCount = 0;
    my $intronCount = 0;
    my $noncodingRNACount = 0;
    my $intergenicCount = 0;
    my $otherCount = 0;
    
    ######## count annotation information
    my $radialSVMCount = 0;
    my $funseqCount = 0;
    my $cnvCount = 0;
    
    while(<RADIAL>){
        chomp;
        my @line;
        my $key;
        @line = split(/\t/, $_);
        $key = "$line[2],$line[3],$line[4],$line[5],$line[6]";
        $radialSVM{$key} = $line[1];
    }
    while(<FUNSEQ>){
        chomp;
        my @line;
        my $key;
        @line = split(/\t/, $_);
        $key = "$line[2],$line[3],$line[4],$line[5],$line[6]";
        $funseq{$key} = $line[1];
    }

    # cnv cannot be processed using key and value
    # create a file to temporarily store

    open(CNVBED, ">$cnvbed") or die "ERROR: cannot open $cnvbed for write:\n";
    
    while(<CNV>){
        chomp;
        my @line;
        my $key;
        my $score;
        @line = split(/\t/, $_);
        $score = $line[1];
        $score =~ /Score=(.*);/;
        $score = $1;
	# chr start end score
	print CNVBED "$line[2]\t$line[3]\t$line[4]\t$score\n";
   }
    close CNVBED;
    

    # get intersect
    !system("$bedtools intersect -a $cnvbed -b $genebed -wa > $cnvfinal") or die "ERROR: cannot find intersect using bedtools, please check whether or not you have installed bedtools\n";
    
    open(CNVFINAL, "$cnvfinal") or die "ERROR: cannot open $cnvfinal for read:\n";
    while(<CNVFINAL>){
	chomp;
	my @line = split("\t", $_);
	# note that this key is different for snv
	my $key = "$line[0],$line[1],$line[2]";             
	$cnv{$key} = $line[3]; 
    }

    while(<EXON>){
        chomp;
        my (@line, @syntax, @content);
        my ($key, $mut, $pro);
        @line = split(/\t/, $_);
        @syntax = split(",", $line[2]);
        @content = split(":", $syntax[0]);
        $mut = $content[3];
        $pro = $content[4];
        $key = "$line[3],$line[4],$line[5],$line[6],$line[7]";
        $exon{$key}{"mutationSyntax"} = $mut;
        $exon{$key}{"proteinSyntax"} = $pro;
        if($line[4] == $line[5]){
            $pointcoding{$key} = 1;
        }
    }
    while(<GENE>){
        chomp;
        my @line;
        my ($key, $gene);                                                            #hash key used for fetch radial SVM score from %radialSVM: mutation->radialSVM
        my ($category, $mutationSyntax, $proteinSyntax, $scoreCategory, $score);
        @line = split(/\t/, $_);
        $gene = $line[1];
        next unless defined $gene;
        $key = "$line[2],$line[3],$line[4],$line[5],$line[6]";
	# note that structural variation key is different !!! chr,start,end;
	my $structKey = "$line[2],$line[3],$line[4]";
        my %cnvScore ; # we also need this hash to store cnv score for each gene 
	
	next unless defined $key;
        next unless defined $structKey;
	

        ####### process gene for noncoding variants
        my @printGene;
        if($gene =~ /(.*?)\(dist=(.*?)\),(.*?)\(dist=(.*?)\)/){
            my $gene1 = $1;
            my $gene2 = $3;
            my $dist1 = $2;
            my $dist2 = $4;
	    if($dist1 eq "NONE" and $dist2 eq "NONE"){
		$printGene[0] = $gene1;
	    }elsif($dist1 eq "NONE"){
		$printGene[0] = $gene2;
	    }elsif($dist2 eq "NONE"){
		$printGene[0] = $gene1;
	    }elsif($dist1 <= $dist2){
                $printGene[0] = $gene1;
            }else{
                $printGene[0] = $gene2;
            }
        }elsif($gene =~ /([A-Z|0-9|-]+?)\(.*\),([A-Z|0-9|-]+?)\(.*\)$/){
            $printGene[0] = $1;
            $printGene[1] = $3;
        }elsif($gene =~ /([A-Z|0-9|-]+?)\(.*\);([A-Z|0-9|-]+?)\(.*\)$/){
            $printGene[0] = $1;
            $printGene[1] = $3;
        }elsif($gene =~ /([A-Z|0-9|-]+?)\(.*\)$/){
            $printGene[0] = $1;
        }elsif($gene =~ /;/ or $gene =~ /,/){
            my @gene = split(/;|,/, $gene);
            for(0..$#gene){
                $printGene[$_] = $gene[$_];
            }
        }else{
            $printGene[0] = $gene;
        }

        
        if($line[0] =~ /^exonic/ || $line[0] =~ /^splicing/ ){
            $exonCount ++;
        }elsif($line[0] =~ /^intron/){
            $intronCount ++;
        }elsif($line[0] =~ /^ncRNA/){
            $noncodingRNACount ++;
        }elsif($line[0] =~ /^intergenic/){
            $intergenicCount ++;
        }else{
            $otherCount ++;
        }
       
        if ($line[3] == $line[4]){
            if(exists $pointcoding{$key}){
                $radialSVMCount ++;
                $category = "point coding";
                if(defined $exon{$key}{"mutationSyntax"}){
                    $mutationSyntax = $exon{$key}{"mutationSyntax"};
                }else{
                    $mutationSyntax = "NA";
                }
                if(defined $exon{$key}{"proteinSyntax"}){
                    $proteinSyntax = $exon{$key}{"proteinSyntax"};
                }else{
                     $proteinSyntax = "NA";
                }
                $scoreCategory = "radial SVM";
                if(exists $radialSVM{$key}){
                    $score = $radialSVM{$key}
                }else{
                    $score = "NA";
                }
            }else{
                $funseqCount ++;
                $category = "point noncoding";
                $mutationSyntax = "NA";
                $proteinSyntax = "NA";
                $scoreCategory = "FunSeq2";
                if(exists $funseq{$key}){
                    $score = $funseq{$key}
                }else{
                    $score = "NA";
                }
            }
        }else{
            $cnvCount ++ ;
            $category = "structural variation";
            if(exists $exon{$key}{"mutationSyntax"} and exists $exon{$key}{"proteinSyntax"}){
                $mutationSyntax = $exon{$key}{"mutationSyntax"};
                $proteinSyntax = $exon{$key}{"proteinSyntax"};
            }else{
                $mutationSyntax = "NA";
                $proteinSyntax = "NA";
            }
            $scoreCategory = "CNV normalized signal";
	    for(0..$#printGene){
		if(exists $cnv{$structKey} and (exists $onc{$printGene[$_]} or exists $sup{$printGene[$_]})){
		    $score = $cnv{$structKey};
		}else{
		    $score = "NA";
		}
		$cnvScore{$printGene[$_]} = $score;
	    }
        }
	
        for(0..$#printGene){
            $icagesMutations{$printGene[$_]}{$key}{"category"} = $category;
            $icagesMutations{$printGene[$_]}{$key}{"mutationSyntax"} = $mutationSyntax;
            $icagesMutations{$printGene[$_]}{$key}{"proteinSyntax"} = $proteinSyntax;
            $icagesMutations{$printGene[$_]}{$key}{"scoreCategory"} = $scoreCategory;
	    if($category eq "structural variation"){
		$icagesMutations{$printGene[$_]}{$key}{"score"} = $cnvScore{$printGene[$_]};
	    }else{
		$icagesMutations{$printGene[$_]}{$key}{"score"} = $score;
	    }
        }
    }
    print OUT "geneName,chrmosomeNumber,start,end,reference,alternative,category,mutationSyntax,proteinSyntax,scoreCategory,mutationScore\n";
    foreach my $gene (sort keys %icagesMutations){
        foreach my $mutation (sort keys %{$icagesMutations{$gene}}){
           if (defined $gene and defined $mutation and defined $icagesMutations{$gene}{$mutation}{"category"} and defined $icagesMutations{$gene}{$mutation}{"mutationSyntax"} and defined $icagesMutations{$gene}{$mutation}{"proteinSyntax"} and defined $icagesMutations{$gene}{$mutation}{"scoreCategory"} and defined $icagesMutations{$gene}{$mutation}{"score"}){
               print OUT "$gene,$mutation,$icagesMutations{$gene}{$mutation}{\"category\"},$icagesMutations{$gene}{$mutation}{\"mutationSyntax\"},$icagesMutations{$gene}{$mutation}{\"proteinSyntax\"},$icagesMutations{$gene}{$mutation}{\"scoreCategory\"},$icagesMutations{$gene}{$mutation}{\"score\"}\n" ;
            }else{
                print "$gene,$mutation,$icagesMutations{$gene}{$mutation}{\"category\"},$icagesMutations{$gene}{$mutation}{\"mutationSyntax\"},$icagesMutations{$gene}{$mutation}{\"proteinSyntax\"},$icagesMutations{$gene}{$mutation}{\"scoreCategory\"},$icagesMutations{$gene}{$mutation}{\"score\"}\n" ;
            }
        }
    }
    
    my $logFile = $annovarInputFile . ".icages.log";
    open(LOG, ">>$logFile") or die "iCAGES: cannot open file $logFile\n";
    
    print LOG "## location information\n";
    print LOG "exonic/splice: $exonCount\n";
    print LOG "intronic: $intronCount\n";
    print LOG "noncoding RNA: $noncodingRNACount\n";
    print LOG "intergenic: $intergenicCount\n";
    print LOG "other: $otherCount\n\n";
    
    print LOG "## annotation information\n";
    print LOG "point coding variants with radialSVM annotation: $radialSVMCount\n";
    print LOG "point noncoding variants with FunSeq2 annotation: $funseqCount\n";
    print LOG "Indels and SVs with CNV signal annotation: $cnvCount\n\n";
}


sub formatConvert{
    # $rawInputFile, $annovarInputFile, $icagesLocation, $tumor, $germline , $id, $bed, $expression
    print "NOTICE: start input file format checking and converting format if needed\n";
    my ($rawInputFile, $annovarInputFile, $icagesLocation );
    my ( $tumor, $germline , $id, $bed, $prefix , $expression); # parameters for vcf conversion
    my  $callConvertToAnnovar;
    my $callvcftools;
    my $formatCheckFirstLine;
    my $isbedFormat = 0; # check whether or not this file is in bed format
    $rawInputFile = shift;
    $annovarInputFile = shift;
    $icagesLocation = shift;
    $tumor = shift;
    $germline = shift;
    $id = shift;
    $bed = shift;
    $expression = shift;
    open(IN, "$rawInputFile") or die "ERROR: cannot open $rawInputFile\n";
    $formatCheckFirstLine = <IN>;
    chomp $formatCheckFirstLine;
    my $multipleSampleCheck = 0;    
    if($formatCheckFirstLine =~ /^#/){
	# check the whole file to see if this is multple sample

	while(<IN>){
	    chomp;
	    my $line = $_;
	    my @line = split;
	    if($line[0] =~ /^#CHROM/){
		if($#line > 9){
		    $multipleSampleCheck = 1;
		}
		last;
	    }
	}
    }else{
	my @line = split(/\t| /, $formatCheckFirstLine);
	if($#line == 2 ){
	    $isbedFormat  = 1;
	}elsif($#line != 2 and defined  $line[3]  and defined $line[4] ){
	    if($line[3] !~ /[a|t|c|g|A|T|C|G|-]+/ or $line[4] !~ /[a|t|c|g|A|T|C|G|-]+/){
		$isbedFormat  = 1;
	    }
	}
    }
    close IN;
    $callConvertToAnnovar = $icagesLocation . "bin/annovar/convert2annovar.pl";
    $callvcftools = $icagesLocation . "bin/vcftools/bin/vcftools";
    if($formatCheckFirstLine =~ /^##fileformat=VCF/){             #VCF
	if($multipleSampleCheck and $tumor eq "NA" and $germline eq "NA" and $id eq "NA"){
	    die "ERROR: your vcf file contains multiple samples please specify a valid sample identifier \n";
	}
	if($tumor ne "NA" and $germline ne "NA"){
	    !system("$callvcftools --recode --vcf $rawInputFile --indv $tumor --out $rawInputFile.$tumor") or die "ERROR: please specify a valid sample identifier for tumor sample\n";
	    !system("$callvcftools --recode --vcf $rawInputFile --indv $germline --out $rawInputFile.$germline") or die "ERROR: please specify a valid sample identifier for germline variants\n";
	    !system("$callConvertToAnnovar -format vcf4 $rawInputFile.$tumor.recode.vcf > $rawInputFile.$tumor.ann") or die "ERROR: cannot execute convert2annovar.pl for converting VCF file\n";
	    !system("$callConvertToAnnovar -format vcf4 $rawInputFile.$germline.recode.vcf > $rawInputFile.$germline.ann") or die "ERROR: cannot execute convert2annovar.pl for converting VCF file\n";
	    !system("cat $rawInputFile.$tumor.ann $rawInputFile.$germline.ann | uniq -u > $annovarInputFile") or die "ERROR: cannot generate somatic variants input file for iCAGES, please double check teh format of your input files\n";
	}elsif($id ne "NA"){
	    !system("$callvcftools --recode --vcf $rawInputFile --indv $id --out $rawInputFile.$id") or die "ERROR: please specify a valid sample identifier for somatic variants of your interest\n";
	    !system("$callConvertToAnnovar -format vcf4 $rawInputFile.$id.recode.vcf > $annovarInputFile") or die "ERROR: cannot execute convert2annovar.pl for converting VCF file\n";
	}else{
	    !system("$callConvertToAnnovar -format vcf4 $rawInputFile > $annovarInputFile") or die "ERROR: cannot execute convert2annovar.pl for converting VCF file\n";
	}
        
    }elsif($isbedFormat){
	# BED
	print "iCAGES: your input file is likely to be a bed file and iCAGES is converting it to ANNOVAR input format\n";
#	!system("$callConvertToAnnovar -format bed $rawInputFile >  $annovarInputFile") or die "ERROR: cannot convert your BED file input into ANNOVAR input format, please double check your input file\n";
	my @cnv;
	open(INPUTCNV, "$rawInputFile") or die;
	my $checkFiveFields = 0;
	while(<INPUTCNV>){
	    chomp;
	    push @cnv, $_;
	    my @line  = split;
	    $checkFiveFields = 1 if $#line == 4;
	}
	close INPUTCNV;
	open(OUTPUTCNV, ">$annovarInputFile") or die;
	for(0..$#cnv){
	    if( $checkFiveFields == 1){
		print OUTPUTCNV "$cnv[$_]\n";
	    }else{
		print OUTPUTCNV "$cnv[$_]\t0\t0\n";
	    }
	}
    }else{                    
	#ANNOVAR
        !system("cp $rawInputFile $annovarInputFile") or die "ERROR: cannot use input file $rawInputFile\n";
    }
    if($bed ne "NA"){
# there is a bug in annovar convert2annovar.pl for bed
#	!system("$callConvertToAnnovar -format bed $bed >  $bed.out") or die "ERROR: cannot convert your BED file input into ANNOVAR input format, please double check your input file\n";
#	!system('awk \'{print $1 "\t" $2 "\t" $3 "\t0\t0"  }\' $bed >  $bed.out') or die "ERROR: cannot convert your BED file input into ANNOVAR input format\n";
	open(ANNBED,  ">>", $annovarInputFile) or die "iCAGES: cannot find converted input file for iCAGES in ANNOVAR input format\n";
	open(BED, "$bed") or die "iCAGES: cannot open ANNOVAR input file generated by your input BED file\n";
	my @bed;
	while(<BED>){
	    chomp;
	    push @bed, $_;
	}
	close BED;
        for(0..$#bed){
	    print ANNBED "$bed[$_]\t0\t0\n";
	}
	close ANNBED;
    }
    if($expression ne "NA"){
	!system("$callConvertToAnnovar -format bed $expression >  $expression.out") or die "ERROR: cannot convert your BED file input into ANNOVAR input format, please double check your input file\n";
        open(ANNBED,  ">>", $annovarInputFile) or die "iCAGES: cannot find converted input file for iCAGES in ANNOVAR input format\n";
        open(EXP, "$expression.out") or die "iCAGES: cannot open ANNOVAR input file generated by your input BED file\n";
        my @bed;
        while(<EXP>){
            chomp;
            push @bed, $_;
        }
	close EXP;
        for(0..$#bed){
            print ANNBED "$bed[$_]\n";
        }
        close ANNBED;
    }
}


sub divideMutation{
    print "NOTICE: start dividing mutations to SNP and structural variation\n";
    my ($annovarInputFile, $snpFile, $cnvFile);
    $annovarInputFile = shift;
    $snpFile = $annovarInputFile . ".snp";
    $cnvFile = $annovarInputFile . ".cnv";
    
    open(OUT, "$annovarInputFile") or die "iCAGES: cannot open input file $annovarInputFile\n";
    open(SNP, ">$snpFile") or die "iCAGES: cannot open input file $snpFile\n";
    open(CNV, ">$cnvFile") or die "iCAGES: cannot open input file $cnvFile\n";
    
    #### add variant information into log file
    my $logFile = $annovarInputFile . ".icages.log";
    open(LOG, ">$logFile") or die "iCAGES: cannot open file $logFile\n";
    my $variantCount = 0;
    my $snvCount = 0;
    my $cnvCount = 0;
    
    while(<OUT>){
        chomp;
        $variantCount ++;
        my $printLine = $_;
        my @line = split(/\t/, $_);
        if ($line[1] == $line[2] and $line[3] ne "-" and $line[4] ne "-"){
            $snvCount ++;
            print SNP "$printLine\n";
        }else{
            $cnvCount ++;
            print CNV "$printLine\n";
        }
    }
    close OUT;
    close SNP;
    close CNV;
    
    print LOG "########### iCAGES Variant Summary ###########\n";
    print LOG "## basic information\n";
    print LOG "Total: $variantCount\n";
    print LOG "SNVs: $snvCount\n";
    print LOG "Indels and Structural variants: $cnvCount\n\n";

}

sub loadDatabase{
    print "NOTICE: start loading databases\n";
    my $DBLocation = shift;
    my $supLocation = $DBLocation . "suppressor.gene";
    my $oncLocation = $DBLocation . "oncogene.gene";
    print "NOTICE: start extracting suppressor genes\n";
    open(SUP, "$supLocation") or die "cannot open $supLocation\n";
    while(<SUP>){
        chomp;
        $sup{$_} = 1;
    }
    close SUP;
    print "NOTICE: start extracting oncogenes\n";
    open(ONC, "$oncLocation") or die "cannot open $oncLocation\n";
    while(<ONC>){
        chomp;
        $onc{$_} = 1;
    }
    close ONC;
}



sub annotateMutation{
    my ($DBLocation, $icagesLocation, $callAnnovar, $annovarInputFile, $snpFile, $cnvFile, $hg);
    $icagesLocation = shift;
    $annovarInputFile = shift;
    $hg = shift;
    $DBLocation = $icagesLocation . "db/";
    $snpFile = $annovarInputFile . ".snp";
    $cnvFile = $annovarInputFile . ".cnv";
    $callAnnovar = $icagesLocation . "bin/annovar/annotate_variation.pl";
    my @children_pids;
    $children_pids[0] = fork();
    if($children_pids[0] == 0){
        print "NOTICE: start to run ANNOVAR region annotation to annotate structural variations or variants associated with CNV changes\n";
	if(-s $cnvFile){
	    !system("$callAnnovar -regionanno -build $hg -out $cnvFile -dbtype cnv $cnvFile $DBLocation -scorecolumn 4 --colsWanted 0") or die "ERROR: cannot call structural varation\n";
	}else{
	    print "NOTICE: CNV file has 0 size\n";
	}
        exit 0;
    }
    $children_pids[1] = fork();
    if($children_pids[1] == 0){
        print "NOTICE: start to run ANNOVAR index function to fetch radial SVM score for each mutation \n";
        if(-s $snpFile){
	    !system("$callAnnovar -filter -out $snpFile -build $hg -dbtype iCAGES $snpFile $DBLocation") or die "ERROR: cannot call icages\n";
        }else{
	    print "NOTICE: SNV file has 0 size\n";
	}
	exit 0;
    }
    $children_pids[2] = fork();
    if($children_pids[2] == 0){
        print "NOTICE: start to run ANNOVAR index function to fetch funseq score for each mutation \n";
	if(-s $snpFile){
	    !system("$callAnnovar -filter -out $snpFile -build $hg -dbtype funseq2 $snpFile $DBLocation") or die "ERROR: cannot call funseq2\n";
	}else{
	    print "NOTICE: SNV file has 0 size, skip funseq score annotation\n";
	}
        exit 0;
    }
    $children_pids[3] = fork();
    if($children_pids[3] == 0){
        print "NOTICE: start annotating each mutaiton using ANNOVAR\n";
        !system("$callAnnovar -out $annovarInputFile -build $hg $annovarInputFile $DBLocation") or die "ERROR: cannot call annovar\n";
        exit 0;
    }
    for (0.. $#children_pids){
        waitpid($children_pids[$_], 0);
    }
}



