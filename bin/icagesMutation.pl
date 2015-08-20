#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);	
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my ($annovarInputFile,$rawInputFile, $inputDir ,$icagesLocation , $tumor ,$germline, $id ,$prefix, $bed, $hg);
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
$nowString = localtime();
$annovarInputFile = &runAnnovar($rawInputFile, $inputDir ,$icagesLocation ,$tumor ,$germline ,$id, $prefix, $bed, $hg );
&processAnnovar($annovarInputFile, $hg);

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
    $rawInputFile = shift;
    $inputDir = shift;
    $icagesLocation = shift;
    $tumor = shift;
    $germline = shift;
    $id = shift;
    $prefix = shift;
    $bed = shift;
    $hg = shift;
    $callAnnotateVariation = $icagesLocation . "bin/annovar/annotate_variation.pl";
    $annovarInputFile = $inputDir . "/" . $prefix . ".annovar";
    $DBLocation = $icagesLocation . "db/";
    &formatConvert($rawInputFile, $annovarInputFile, $icagesLocation, $tumor, $germline , $id , $bed );
    &divideMutation($annovarInputFile);
    &loadDatabase($DBLocation);
    &annotateMutation($icagesLocation, $annovarInputFile, $hg);
    return $annovarInputFile;
}



sub processAnnovar{
    print "NOTICE: start processing output from ANNOVAR\n";
    my $annovarInputFile = shift;
    my $hg = shift;
    my $annovarVariantFunction = $annovarInputFile . ".variant_function";
    my $annovarExonVariantFunction = $annovarInputFile . ".exonic_variant_function";
    my $annovarRadialSVM = $annovarInputFile . ".snp." . $hg . "_iCAGES_dropped";
    my $annovarCNV = $annovarInputFile . ".cnv." . $hg . "_cnv";
    my $annovarFunseq2 = $annovarInputFile . ".snp." . $hg . "_funseq2_dropped";
    my $icagesMutations = $annovarInputFile . ".icagesMutations.csv";
    open(GENE, "$annovarVariantFunction") or die "ERROR: cannot open file $annovarVariantFunction\n";
    open(EXON, "$annovarExonVariantFunction") or die "ERROR: cannot open file $annovarExonVariantFunction\n";
    open(CNV, "$annovarCNV") or die "ERROR: cannot open file $annovarCNV\n";
    open(RADIAL, "$annovarRadialSVM") or die "ERROR: cannot open file $annovarRadialSVM\n";
    open(FUNSEQ, "$annovarFunseq2") or die "ERROR: cannot open file $annovarFunseq2\n";
    open(OUT, ">$icagesMutations") or die "ERROR: cannot open file $icagesMutations\n";
    my (%radialSVM, %funseq, %cnv, %exon);
    my (%pointcoding);
    my %icagesMutations;
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
    while(<CNV>){
        chomp;
        my @line;
        my $key;
        my $score;
        @line = split(/\t/, $_);
        $score = $line[1];
        $score =~ /Score=(.*);/;
        $score = $1;
        $key = "$line[2],$line[3],$line[4],$line[5],$line[6]";
        $cnv{$key} = $score;
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
        next unless defined $key;
        if ($line[3] == $line[4]){
            if(exists $pointcoding{$key}){
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
            $category = "structural variation";
            if(exists $exon{$key}{"mutationSyntax"} and exists $exon{$key}{"proteinSyntax"}){
                $mutationSyntax = $exon{$key}{"mutationSyntax"};
                $proteinSyntax = $exon{$key}{"proteinSyntax"};
            }else{
                $mutationSyntax = "NA";
                $proteinSyntax = "NA";
            }
            $scoreCategory = "CNV normalized signal";
            if(exists $cnv{$key} and exists $sup{$gene}){
                $score = $cnv{$key}
            }else{
                $score = "NA";
            }
        }
        $icagesMutations{$gene}{$key}{"category"} = $category;
        $icagesMutations{$gene}{$key}{"mutationSyntax"} = $mutationSyntax;
        $icagesMutations{$gene}{$key}{"proteinSyntax"} = $proteinSyntax;
        $icagesMutations{$gene}{$key}{"scoreCategory"} = $scoreCategory;
        $icagesMutations{$gene}{$key}{"score"} = $score;
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
}


sub formatConvert{
    # $rawInputFile, $annovarInputFile, $icagesLocation, $tumor, $germline , $id, $bed
    print "NOTICE: start input file format checking and converting format if needed\n";
    my ($rawInputFile, $annovarInputFile);
    my ( $tumor, $germline , $id, $prefix ); # parameters for vcf conversion
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
    open(IN, "$rawInputFile") or die "ERROR: cannot open $rawInputFile\n";
    $formatCheckFirstLine = <IN>;
    chomp $formatCheckFirstLine;
    my $multipleSampleCheck = 0;
    while(<IN>){
	chomp;
	my $line = $_;
        my @line = split;
	if($line[0] =~ /^#CHROM/){
	    if($#line > 9){
		$multipleSampleCheck = 1;
	    }
	    last;
	}elsif($#line == 2 ){
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
	!system("$callConvertToAnnovar -format bed $rawInputFile >  $annovarInputFile") or die "ERROR: cannot convert your BED file input into ANNOVAR input format, please double check your input file\n";
    }else{                    
	#ANNOVAR
        !system("cp $rawInputFile $annovarInputFile") or die "ERROR: cannot use input file $rawInputFile\n";
    }
    if($bed ne "NA"){
	!system("$callConvertToAnnovar -format bed $bed >  $bed.out") or die "ERROR: cannot convert your BED file input into ANNOVAR input format, please double check your input file\n";
	open(ANNBED,  ">>", $annovarInputFile) or die "iCAGES: cannot find converted input file for iCAGES in ANNOVAR input format\n";
	open(BED, "$bed.out") or die "iCAGES: cannot open ANNOVAR input file generated by your input BED file\n";
	my @bed;
	while(<BED>){
	    chomp;
	    push @bed, $_;
	}
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
    while(<OUT>){
        chomp;
        my $printLine = $_;
        my @line = split(/\t/, $_);
        if ($line[1] == $line[2] and $line[3] ne "-" and $line[4] ne "-"){
            print SNP "$printLine\n";
        }else{
            print CNV "$printLine\n";
        }
    }
    close OUT;
    close SNP;
    close CNV;

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
        print "NOTICE: start to run ANNOVAR region annotation to annotate structural variations or variants associated with LOF changes\n";
        !system("$callAnnovar -regionanno -build $hg -out $cnvFile -dbtype cnv $cnvFile $DBLocation -scorecolumn 4 --colsWanted 0") or die "ERROR: cannot call structural varation\n";
        exit 0;
    }
    $children_pids[1] = fork();
    if($children_pids[1] == 0){
        print "NOTICE: start to run ANNOVAR index function to fetch radial SVM score for each mutation \n";
        !system("$callAnnovar -filter -out $snpFile -build $hg -dbtype iCAGES $snpFile $DBLocation") or die "ERROR: cannot call icages\n";
        exit 0;
    }
    $children_pids[2] = fork();
    if($children_pids[2] == 0){
        print "NOTICE: start to run ANNOVAR index function to fetch funseq score for each mutation \n";
        !system("$callAnnovar -filter -out $snpFile -build $hg -dbtype funseq2 $snpFile $DBLocation") or die "ERROR: cannot call funseq2\n";
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



