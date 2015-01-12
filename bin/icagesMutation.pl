#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);	
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my ($annovarInputFile, $rawInputFile, $icagesLocation);
my (%sup, %onc);

######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

$rawInputFile = $ARGV[0];
$icagesLocation = $ARGV[1];
$annovarInputFile = &runAnnovar($rawInputFile, $icagesLocation);
&processAnnovar($annovarInputFile, $rawInputFile);

######################################################################################################################################
############################################################# subroutines ############################################################
######################################################################################################################################

sub runAnnovar {
    print "NOTICE: start runing iCAGES packge at $nowString\n";
    my ($rawInputFile, $annovarInputFile);                                                                      #ANNOVAR input files
    my ($icagesLocation, $callAnnotateVariation);                                                               #ANNOVAR commands
    my ($radialSVMDB, $radialSVMIndex, $funseq2DB, $funseq2Index, $refGeneDB, $refGeneIndex, $cnvDB);           #ANNOVAR DB files: iCAGES score (index), refGene (fasta), dbSNP
    my ($log, $annovarLog);
    my $nowString = localtime;                                                                                  #SYSTEM local time
    $rawInputFile = shift;
    $icagesLocation = shift;
    $callAnnotateVariation = $icagesLocation . "bin/annovar/annotate_variation.pl";
    $annovarInputFile = $rawInputFile . ".annovar";
    $DBLocation = $icagesLocation . "db/";
    &formatConvert($rawInputFile, $annovarInputFile, $icagesLocation);
    &divideMutation($annovarInputFile);
    &loadDatabase($DBLocation);
    &annotateMutation($DBLocation, $annovarInputFile);
    return $annovarInputFile;
}



sub processAnnovar{
    print "NOTICE: start processing output from ANNOVAR\n";
    my $annovarInputFile = shift;
    my $rawInputFile = shift;
    my $annovarVariantFunction = $annovarInputFile . ".variant_function";
    my $annovarExonVariantFunction = $annovarInputFile . ".exonic_variant_function";
    my $annovarRadialSVM = $annovarInputFile . ".snp.hg19_iCAGES_dropped";
    my $annovarCNV = $annovarInputFile . ".cnv.hg19_cnv";
    my $annovarFunseq2 = $annovarInputFile . ".snp.hg19_funseq2_dropped";
    my $icagesMutations = $rawInputFile . ".icagesMutations.csv";
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
        $key = "$line[3]:$line[4]:$line[5]:$line[6]:$line[7]";
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
        $key = "$line[2],$line[3],$line[4],$line[5],$line[6]";
        $gene = $line[1];
        next unless defined $gene;
        if ($line[5] == $line[6]){
            if(exists $pointcoding{$key}){
                $category = "point coding";
                $mutationSyntax = $exon{$key}{"mutationSyntax"};
                $proteinSyntax = $exon{$key}{"proteinSyntax"};
                $scoreCategory = "radial SVM";
                if(exists $radialSVM{$key}){
                    $score = $radialSVM{$key}
                }else{
                    $score = "NA";
                }
            }else{
                $category = "point coding";
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
            if(exists $exon{$key}){
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
    
    foreach my $gene (sort keys %icagesMutations){
        foreach my $mutation (sort keys %{$icagesMutation{$gene}}){
            print OUT "$gene,$key,$icagesMutations{$gene}{$key}{\"category\"},$icagesMutations{$gene}{$key}{\"mutationSyntax\"},$icagesMutations{$gene}{$key}{\"proteinSyntax\"},$icagesMutations{$gene}{$key}{\"scoreCategory\"},$icagesMutations{$gene}{$key}{\"score\"}\n";
        }
    }
}


sub formatConvert{
    print "NOTICE: start input file format checking and converting format if needed\n";
    my ($rawInputFile, $annovarInputFile);
    my  $callConvertToAnnovar;
    my $formatCheckFirstLine;
    $rawInputFile = shift;
    $annovarInputFile = shift;
    $icagesLocation = shift;
    open(IN, "$rawInputFile") or die "ERROR: cannot open $rawInputFile\n";
    $formatCheckFirstLine = <IN>;
    chomp $formatCheckFirstLine;
    close IN;
    $callConvertToAnnovar = $icagesLocation . "bin/DGIdb/convert2annovar.pl";
    if($formatCheckFirstLine =~ /^##fileformat=VCF/){                                                               #VCF
        !system("$callConvertToAnnovar -format vcf4 $rawInputFile > $annovarInputFile") or die "ERROR: cannot execute convert2annovar.pl for converting VCF file\n";
    }else{                                                                                                          #ANNOVAR
        !system("cp $rawInputFile $annovarInputFile") or die "ERROR: cannot use input file $rawInputFile\n";
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
    my ($DBLocation, $annovarInputFile, $snpFile, $cnvFile);
    $DBLocation = shift;
    $annovarInputFile = shift;
    $snpFile = $annovarInputFile . ".snp";
    $cnvFile = $annovarInputFile . ".cnv";
    my @children_pids;
    $children_pids[0] = fork();
    if($children_pids[0] == 0){
        print "NOTICE: start to run ANNOVAR region annotation to annotate structural variations or variants associated with LOF changes\n";
        push (@log, "iCAGES: start to run ANNOVAR region annotation to annotate structural variations or variants associated with LOF changes");
        !system("$callAnnovar -regionanno -build hg19 -out $cnvFile -dbtype cnv $cnvFile $DBLocation -scorecolumn 4 --colsWanted 0") or die "ERROR: cannot call structural varation\n";
        exit 0;
    }
    $children_pids[1] = fork();
    if($children_pids[1] == 0){
        print "NOTICE: start to run ANNOVAR index function to fetch radial SVM score for each mutation \n";
        push (@log, "iCAGES: start to run ANNOVAR index function to fetch radial SVM score for each mutation");
        !system("$callAnnovar -filter -out $snpFile -build hg19 -dbtype iCAGES $snpFile $DBLocation") or die "ERROR: cannot call icages\n";
        exit 0;
    }
    $children_pids[2] = fork();
    if($children_pids[2] == 0){
        print "NOTICE: start to run ANNOVAR index function to fetch funseq score for each mutation \n";
        push (@log, "iCAGES: start to run ANNOVAR index function to fetch funseq score for each mutation");
        !system("$callAnnovar -filter -out $snpFile -build hg19 -dbtype funseq2 $snpFile $DBLocation") or die "ERROR: cannot call funseq2\n";
        exit 0;
    }
    $children_pids[3] = fork();
    if($children_pids[3] == 0){
        print "NOTICE: start annotating each mutaiton using ANNOVAR\n";
        push (@log, "iCAGES: start annotating each mutaiton using ANNOVAR");
        !system("$callAnnovar -out $outputFile -build hg19 $outputFile $DBLocation") or die "ERROR: cannot call annovar\n";
        exit 0;
    }
    for (0.. $#children_pids){
        waitpid($children_pids[$_], 0);
    }
}



