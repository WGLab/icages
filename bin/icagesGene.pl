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
	next unless defined $line[0];
        if(defined $icagesGenes{$line[0]}{$line[9]}){
            if($line[10] eq "NA" or  $icagesGenes{$line[0]}{$line[9]} eq "NA"){
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

	next if $cnv==0 and $funseq==0 and $radialSVM==0;
	
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
            $icagesGene = 1.36771 -0.42382 * $radialSVM -0.93636 * $cnv  -0.15453 * $funseq -2.43656 * $phenolyzer;
        }elsif($subtype eq "BLCA"){
            $icagesGene = 1.058757 -0.038520 * $radialSVM -0.398104 * $cnv -0.003739 * $funseq -0.909022 * $phenolyzer;
        }elsif($subtype eq "BRCA"){
            $icagesGene = 1.093442 -0.075614 * $radialSVM  -0.479842* $cnv -0.010534 * $funseq -1.022521 * $phenolyzer;
        }elsif($subtype eq "CESC"){
            $icagesGene = 1.16180 -0.14885 * $radialSVM -1.32652 * $cnv -0.05598 * $funseq -1.32015 * $phenolyzer;
        }elsif($subtype eq "CHOL"){
            $icagesGene = -0.14641 + 0.18845 * $radialSVM + 0.50381 * $cnv +0.13216 * $funseq + 0.61047 * $phenolyzer;
        }elsif($subtype eq "COADREAD") {
	    $icagesGene = -0.062618 + 0.041399 * $radialSVM + 1.315360 * $cnv + 0.008455 * $funseq + 0.990089 * $phenolyzer;
	}elsif($subtype eq "COAD"){
            $icagesGene = -0.091516 + 0.069962 * $radialSVM + 1.101200 * $cnv  + 0.030804 * $funseq + 1.044726 * $phenolyzer;
        }elsif($subtype eq "DLBC") {
	    $icagesGene = -0.12036 + 0.13884 * $radialSVM + 0.53773 * $cnv  + 0.01512 * $funseq + 0.85130 * $phenolyzer;
	}elsif($subtype eq "ESCA"){
            $icagesGene = -0.13015 + 0.11817 * $radialSVM + 0.82129 * $cnv + 0.01265 * $funseq + 1.07817 * $phenolyzer;
        }elsif($subtype eq "GBM"){
            $icagesGene = -0.12775 + 0.12841 * $radialSVM + 0.45996 * $cnv + 0.02754 * $funseq + 0.92774 * $phenolyzer;
        }elsif($subtype eq "GBMLGG"){
	    $icagesGene = -0.07477 + 0.05322 * $radialSVM + 1.18877 * $cnv -0.01261 * $funseq + 1.08833 * $phenolyzer;
	} elsif($subtype eq "HNSC"){
            $icagesGene = 1.058015 -0.037751 * $radialSVM -0.581648 * $cnv -0.004731 * $funseq -0.946285 * $phenolyzer;
        }elsif($subtype eq "KICH"){
            $icagesGene = -0.12298 + 0.16503 * $radialSVM + 0.29233 * $cnv + 0.10256 * $funseq + 0.46457 * $phenolyzer;
        }elsif($subtype eq "KIPAN") {
	    $icagesGene = 1.093222 -0.073644 * $radialSVM -0.345129 * $cnv - 0.028881 * $funseq -1.146782 * $phenolyzer;
	}elsif($subtype eq "KIRC"){
            $icagesGene = 1.13778 -0.12890 * $radialSVM -0.54436 * $cnv -0.01729 * $funseq -1.44148 * $phenolyzer;
        }elsif($subtype eq "KIRP"){
            $icagesGene = -0.07488 + 0.06861 * $radialSVM + 1.23478 * $cnv +0.01712 * $funseq + 0.96863 * $phenolyzer;
        }elsif($subtype eq "LAML"){
            $icagesGene = -0.05995 + 0.15925 * $radialSVM + 0.12114 * $cnv + 0.02738 * $funseq + 0.02367 * $phenolyzer;
        }elsif($subtype eq "LGG"){
            $icagesGene = -0.11611 + 0.10186 * $radialSVM + 0.95189 * $cnv + 0.04495 * $funseq + 1.11112 * $phenolyzer;
        }elsif($subtype eq "LIHC"){
            $icagesGene = 1.113955 -0.100192 * $radialSVM -0.773660 * $cnv  -0.047867 * $funseq -1.076089 * $phenolyzer;
        }elsif($subtype eq "LUAD") {
	    $icagesGene = -0.096116 + 0.074096 * $radialSVM + 1.462135 * $cnv + 0.032892 * $funseq + 1.035469 * $phenolyzer;
	}elsif($subtype eq "LUSC"){
            $icagesGene = 1.13182 -0.11578 * $radialSVM -0.53788 * $cnv -0.03912 * $funseq -1.15891 * $phenolyzer;
        }elsif($subtype eq "OV"){
            $icagesGene = 1.18324 -0.17481 * $radialSVM -0.72561 * $cnv -0.06436 * $funseq -1.42091 * $phenolyzer;
        }elsif($subtype eq "PAAD"){
            $icagesGene = 1.16658 -0.15211 * $radialSVM -0.87784 * $cnv -0.08511 * $funseq -1.33967 * $phenolyzer;
        }elsif($subtype eq "PCPG"){
            $icagesGene = -0.12564 + 0.19600 * $radialSVM + 0.08097 * $cnv  + 0.06315 * $funseq + 0.28951 * $phenolyzer;
        }elsif($subtype eq "PRAD"){
            $icagesGene = 1.124113 -0.105090 * $radialSVM -0.299424 * $cnv -0.001517 * $funseq -1.375397 * $phenolyzer;
        }elsif($subtype eq "READ") {
	    $icagesGene = -0.12973 + 0.12316 * $radialSVM + 1.52412 * $cnv +  0.03016 * $funseq + 0.98996 * $phenolyzer;
	}elsif($subtype eq "SARC"){
            $icagesGene = -0.077275 + 0.081598 * $radialSVM + 0.562010 * $cnv -0.006907 * $funseq + 0.772565 * $phenolyzer;
        }elsif($subtype eq "SKCM"){
            $icagesGene = 1.076621 -0.051706 * $radialSVM -0.283882 * $cnv -0.015553 * $funseq -1.104187 * $phenolyzer;
        }elsif($subtype eq "STAD"){
            $icagesGene = 1.028806 -0.010729 * $radialSVM -0.753154 * $cnv + 0.007523 * $funseq -0.807513 * $phenolyzer;
        }elsif($subtype eq "STES") {
	    $icagesGene = 1.0288425 -0.0127624 * $radialSVM -0.5998427 * $cnv + 0.0004369 * $funseq -0.7377184 * $phenolyzer;
	}elsif($subtype eq "TGCT"){
            $icagesGene = -0.12173 + 0.14025 * $radialSVM -0.01375 * $cnv + 0.01167 * $funseq + 0.73381 * $phenolyzer;
        }elsif($subtype eq "THCA"){
            $icagesGene = -0.106440 + 0.125886 * $radialSVM + 0.753607 * $cnv + 0.008624 * $funseq + 0.705268 * $phenolyzer;
        }elsif($subtype eq "THYM"){
            $icagesGene = -0.07786 + 0.15530 * $radialSVM + 0.25656 * $cnv + 0.04185 * $funseq + 0.15939 * $phenolyzer;
        }elsif($subtype eq "UCEC"){
            $icagesGene = -0.064307 + 0.038792 * $radialSVM + 0.992520 * $cnv + 0.006978 * $funseq + 1.058414 * $phenolyzer;
        }elsif($subtype eq "UCS"){
            $icagesGene = -0.132534 + 0.153048 * $radialSVM + 0.500817 * $cnv -0.004412 * $funseq + 0.766499 * $phenolyzer;
        }else{
	    $icagesGene = -8.0124 + 5.5577 * $radialSVM + 267.3371 * $cnv + 0.3741 * $funseq + 10.0949 * $phenolyzer;
#            $icagesGene = 0.455825 + 0.053429 * $radialSVM + 0.267669 * $cnv + 0.054989 * $funseq  -0.117 * $phenolyzer;
        }
        $icagesGene = 1/(1+exp(-$icagesGene));
        $icagesPrint{$gene}{"score"} = $icagesGene;
	# add driver information requested by user
	my $driver = "No";
	if ($icagesGene >= 0.11) {
	    $driver = "Yes";
	}
        $icagesPrint{$gene}{"content"} = "$gene,$category,$radialSVM,$funseq,$cnv,$phenolyzer,$icagesGene,$driver";
    };
    print GENES "geneName,category,radialSVM,funseq,cnv,phenolyzer,icagesGeneScore,driver\n";
    foreach my $gene (sort {$icagesPrint{$b}{"score"} <=> $icagesPrint{$a}{"score"}} keys %icagesPrint){
        print GENES "$icagesPrint{$gene}{\"content\"}\n";
    }
    
    my $logFile = $rawInputFile . $prefix  . ".annovar.icages.log";
    open(LOG, ">>$logFile") or die "iCAGES: cannot open file $logFile\n";
    
    print LOG "########### iCAGES Gene Summary ###########\n";
    print LOG "## basic information\n";
    print LOG "Total: $geneCount\n";
    print LOG "Cancer Gene Census Gene: $cgcCount\n";
    print LOG "KEGG Pathway Gene: $keggCount\n\n";
    
}


sub loadDatabase {
    print "NOTICE: start loading CGC, KEGG cancer pathway databases\n";
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



