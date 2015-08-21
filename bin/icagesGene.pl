#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);	
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my ($rawInputFile, $icagesLocation, $subtype, $prefix);
my %phenolyzer;

######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

$rawInputFile = $ARGV[0];
$icagesLocation = $ARGV[1];
$subtype = $ARGV[2];
$prefix = $ARGV[3];
%phenolyzer = &loadPhenolyzer($icagesLocation);
&processMutation($rawInputFile, $icagesLocation, \%phenolyzer, $subtype, $prefix);


######################################################################################################################################
############################################################# subroutines ############################################################
######################################################################################################################################

sub loadPhenolyzer {
    print "NOTICE: start loading Phenolyzer\n";
    my %phenolyzer;
    my ($icagesLocation, $DBLocation, $phenolyzerDB);
    $icagesLocation = shift;
    $DBLocation = $icagesLocation . "db/";
    $phenolyzerDB = $DBLocation . "phenolyzer.score";
    open(PHE, "$phenolyzerDB") or die "ERROR: cannot open $phenolyzerDB\n";
    while(<PHE>){
        chomp;
        my @line = split("\t", $_);
        $phenolyzer{$line[0]} = $line[1];
    }
    return %phenolyzer;
}


sub processMutation{
    print "NOTICE: start process mutation files from iCAGES layer one\n";
    my ($rawInputFile, $icagesLocation, $DBLocation, $icagesMutations, $icagesGenes, $ref, $subtype);
    my (%phenolyzer, %icagesGenes, %icagesPrint);
    my (%cgc, %kegg);
    my ($cgcRef, $keggRef);
    $rawInputFile = shift;
    $icagesLocation = shift;
    $ref = shift;
    $subtype = shift;
    $prefix = shift;
    %phenolyzer = %{$ref};
    $DBLocation = $icagesLocation . "db/";
    ($cgcRef, $keggRef) = loadDatabase($DBLocation);
    %cgc = %{$cgcRef};
    %kegg = %{$keggRef};
    $icagesMutations = $rawInputFile . $prefix . ".annovar.icagesMutations.csv";
    $icagesGenes = $rawInputFile . $prefix . ".annovar.icagesGenes.csv";
    open(MUTATIONS, "$icagesMutations") or die "ERROR: cannot open $icagesMutations\n";
    my $header = <MUTATIONS>;
    open(GENES, ">$icagesGenes") or die "ERROR: cannot open $icagesGenes\n";
    while(<MUTATIONS>){
        chomp;
        my @line = split(",", $_);
        if(defined $icagesGenes{$line[0]}{$line[9]}){
            if($line[10] eq "NA"){
                next;
            }else{
                $icagesGenes{$line[0]}{$line[9]} = max($line[10], $icagesGenes{$line[0]}{$line[9]});
            };
        }else{
            $icagesGenes{$line[0]}{$line[9]} = $line[10];
        };
    };
    
    ####### count genes
    my $geneCount = 0;
    my $cgcCount = 0;
    my $keggCount = 0;
    

    foreach my $gene (sort keys %icagesGenes){
        $geneCount ++;
        my ($radialSVM, $funseq, $cnv, $phenolyzer, $icagesGene, $category);
        if (exists $icagesGenes{$gene}{"radial SVM"} and $icagesGenes{$gene}{"radial SVM"} ne "NA"){
            $radialSVM = $icagesGenes{$gene}{"radial SVM"};
        }else{
            $radialSVM = 0;
        };
        if (exists $icagesGenes{$gene}{"FunSeq2"} and $icagesGenes{$gene}{"FunSeq2"} ne "NA"){
            $funseq = $icagesGenes{$gene}{"FunSeq2"} ;
        }else{
            $funseq = 0;
        };
        if (exists $icagesGenes{$gene}{"CNV normalized signal"} and $icagesGenes{$gene}{"CNV normalized signal"} ne "NA"){
            $cnv = $icagesGenes{$gene}{"CNV normalized signal"} ;
        }else{
            $cnv =  0;
        };
        if (exists $phenolyzer{$gene}){
            $phenolyzer = $phenolyzer{$gene};
        }else{
            $phenolyzer = 0;
        };
        if (exists $cgc{$gene}){
            $cgcCount ++;
            $category = "Cancer Gene Census";
        }elsif(exists $kegg{$gene}){
            $keggCount ++;
            $category = "KEGG Cancer Pathway";
        }else{
            $category = "Other Category";
        }
        if($subtype eq "ACC"){
            $icagesGene = 0.030917  - 0.005262 * $radialSVM + 0 * $cnv  - 0.021463 * $funseq + 0.042777 * $phenolyzer;
        }elsif($subtype eq "BLCA"){
            $icagesGene = 0.0003761 + 0.0137814 * $radialSVM + 0 * $cnv + 0.0009165 * $funseq + 0.1323313 * $phenolyzer;
        }elsif($subtype eq "BRCA"){
            $icagesGene = 0.005993 + 0.021377 * $radialSVM + 18.071896 * $cnv + 0.001938 * $funseq + 0.123393 * $phenolyzer;
        }elsif($subtype eq "CESC"){
            $icagesGene = 0.0021918 + 0.0190018 * $radialSVM + 0 * $cnv + 0.0021887 * $funseq + 0.0730895 * $phenolyzer;
        }elsif($subtype eq "CHOL"){
            $icagesGene = 0.036735 + 0.019162 * $radialSVM + 0 * $cnv - 0.018391 * $funseq + 0.122516 * $phenolyzer;
        }elsif($subtype eq "COAD"){
            $icagesGene = 0.007748 + 0.007186 * $radialSVM + 0 * $cnv  - 0.006699 * $funseq + 0.087172 * $phenolyzer;
        }elsif($subtype eq "ESCA"){
            $icagesGene = 0.003368 + 0.038223 * $radialSVM + 0 * $cnv + 0.003475 * $funseq + 0.114999 * $phenolyzer;
        }elsif($subtype eq "GBM"){
            $icagesGene = 0.009045 + 0.018875 * $radialSVM + 17.961815 * $cnv  -0.005318 * $funseq + 0.125607 * $phenolyzer;
        }elsif($subtype eq "HNSC"){
            $icagesGene = 0.0008902 + 0.0153755 * $radialSVM + 0 * $cnv + 0.0068528 * $funseq + 0.1023613 * $phenolyzer;
        }elsif($subtype eq "KICH"){
            $icagesGene = 0.025087 + 0.010040 * $radialSVM + 0 * $cnv + 0.026005 * $funseq + 0.026283 * $phenolyzer;
        }elsif($subtype eq "KIRC"){
            $icagesGene = 0.040234 + 0.037690 * $radialSVM + 0 * $cnv + 0 * $funseq + 0.153110 * $phenolyzer;
        }elsif($subtype eq "KIRP"){
            $icagesGene = 0.0227005 + 0.0322180 * $radialSVM + 0 * $cnv  - 0.0005714 * $funseq + 0.0730942 * $phenolyzer;
        }elsif($subtype eq "LAML"){
            $icagesGene = 0.003622 + 0.096429 * $radialSVM + 0 * $cnv + 0.012641 * $funseq + 0.023935 * $phenolyzer;
        }elsif($subtype eq "LGG"){
            $icagesGene = 0.010569 + 0.011021 * $radialSVM + 0 * $cnv - 0.002998 * $funseq + 0.093497 * $phenolyzer;
        }elsif($subtype eq "LIHC"){
            $icagesGene = 0.0015797 + 0.0042728 * $radialSVM + 0 * $cnv  - 0.0006700 * $funseq + 0.0498169 * $phenolyzer;
        }elsif($subtype eq "LUSC"){
            $icagesGene = 0.002969 + 0.019309 * $radialSVM + 0 * $cnv  - 0.002526 * $funseq + 0.022571 * $phenolyzer;
        }elsif($subtype eq "OV"){
            $icagesGene = 0.004591 + 0.029594 * $radialSVM + 0 * $cnv -0.003883 * $funseq + 0.146881 * $phenolyzer;
        }elsif($subtype eq "PAAD"){
            $icagesGene = 0.007687 + 0.011616 * $radialSVM + 0 * $cnv - 0.006985 * $funseq + 0.098622 * $phenolyzer;
        }elsif($subtype eq "PCPG"){
            $icagesGene = 0.076434 - 0.021682 * $radialSVM + 0 * $cnv  -0.052805 * $funseq + 0.059079 * $phenolyzer;
        }elsif($subtype eq "PRAD"){
            $icagesGene = 0.022258 + 0.017408 * $radialSVM + 0 * $cnv - 0.020219 * $funseq + 0.071447 * $phenolyzer;
        }elsif($subtype eq "SARC"){
            $icagesGene = 0.004694 + 0.026205 * $radialSVM + 0 * $cnv + 0.005975 * $funseq + 0.086695 * $phenolyzer;
        }elsif($subtype eq "SKCM"){
            $icagesGene = 0.0001238 + 0.0023076 * $radialSVM + 0 * $cnv -0.0003054 * $funseq + 0.0535681 * $phenolyzer;
        }elsif($subtype eq "STAD"){
            $icagesGene = 0.0024410 + 0.0307445 * $radialSVM + 0 * $cnv + 0.0003355 * $funseq + 0.1256780 * $phenolyzer;
        }elsif($subtype eq "TGCT"){
            $icagesGene = 0.011837 + 0.014045 * $radialSVM + 0 * $cnv -0.004939 * $funseq + 0.030369 * $phenolyzer;
        }elsif($subtype eq "THCA"){
            $icagesGene = 0.017888 + 0.009046 * $radialSVM + 0 * $cnv -0.019163 * $funseq + 0.171821 * $phenolyzer;
        }elsif($subtype eq "THYM"){
            $icagesGene = 0.008019 + 0.017614 * $radialSVM + 0 * $cnv - 0.005964 * $funseq + 0.031325 * $phenolyzer;
        }elsif($subtype eq "UCEC"){
            $icagesGene = 0.000786 + 0.021023 * $radialSVM + 18.530030 * $cnv + 0.004978 * $funseq + 0.093840 * $phenolyzer;
        }elsif($subtype eq "UCS"){
            $icagesGene = 0.015620 + 0.005112 * $radialSVM + 0 * $cnv -0.014817 * $funseq + 0.092391 * $phenolyzer;
        }elsif($subtype eq "UVM"){
            $icagesGene = 0.060450 + 0.026238 * $radialSVM + 0 * $cnv - 0.058628 * $funseq + 0.164154 * $phenolyzer;
        }else{
            $icagesGene = -8.0124 + 5.5577 * $radialSVM + 267.3371 * $cnv + 0.3741 * $funseq + 10.0949 * $phenolyzer;
        }
        $icagesGene = 1/(1+exp(-$icagesGene));
        $icagesPrint{$gene}{"score"} = $icagesGene;
        $icagesPrint{$gene}{"content"} = "$gene,$category,$radialSVM,$funseq,$cnv,$phenolyzer,$icagesGene";
    };
    print GENES "geneName,category,radialSVM,funseq,cnv,phenolyzer,icagesGeneScore\n";
    foreach my $gene (sort {$icagesPrint{$b}{"score"} <=> $icagesPrint{$a}{"score"}} keys %icagesPrint){
        print GENES "$icagesPrint{$gene}{\"content\"}\n";
    }
    
    my $logFile = $rawInputFile . $prefix  . ".icages.log";
    open(LOG, ">>$logFile") or die "iCAGES: cannot open file $logFile\n";
    
    print LOG "########### iCAGES Gene Summary ###########\n";
    print LOG "## basic information\n";
    print LOG "Total: $geneCount\n";
    print LOG "Cancer Gene Census Gene: $cgcCount\n";
    print LOG "KEGG Pathway Gene: $keggCount\n\n";
    
}


sub loadDatabase {
    print "NOTICE: start loading CGC, KEGG cancer pathway databases";
    my ($DBLocation, $cgcFile, $keggFile);
    my (%cgc, %kegg);
    $DBLocation = shift;
    $cgcFile = $DBLocation . "cgc.gene";
    $keggFile = $DBLocation . "kegg.gene";
    open(CGC, "$cgcFile") or die ;
    open(KEGG, "$keggFile") or die ;
    while(<CGC>){
        chomp;
        my @line = split;
        $cgc{$line[0]} = 1;
    }
    while(<KEGG>){
        chomp;
        my @line = split;
        $kegg{$line[0]} = 1;
    }
    close CGC;
    close KEGG;
    return (\%cgc, \%kegg);
}



