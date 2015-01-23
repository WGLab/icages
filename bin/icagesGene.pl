#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);	
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my ($rawInputFile, $icagesLocation);
my %phenolyzer;

######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

$rawInputFile = $ARGV[0];
$icagesLocation = $ARGV[1];
%phenolyzer = &loadPhenolyzer($icagesLocation);
&processMutation($rawInputFile, $icagesLocation, \%phenolyzer);


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
    my ($rawInputFile, $icagesLocation, $icagesMutations, $icagesGenes, $ref);
    my (%phenolyzer, %icagesGenes, %icagesPrint);
    $rawInputFile = shift;
    $icagesLocation = shift;
    $ref = shift;
    %phenolyzer = %{$ref};
    $icagesMutations = $rawInputFile . ".icagesMutations.csv";
    $icagesGenes = $rawInputFile . ".icagesGenes.csv";
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
    foreach my $gene (sort keys %icagesGenes){
        my ($radialSVM, $funseq, $cnv, $phenolyzer, $icagesGene);
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
        $icagesGene = -8.0124 + 5.5577 * $radialSVM + 267.3371 * $cnv + 0.3741 * $funseq + 10.0949 * $phenolyzer;
        $icagesGene = 1/(1+exp(-$icagesGene));
        $icagesPrint{$gene}{"score"} = $icagesGene;
        $icagesPrint{$gene}{"content"} = "$gene,$radialSVM,$funseq,$cnv,$phenolyzer,$icagesGene";
    };
    print GENES "gene,radialSVM,funseq,cnv,phenolyzer,icagesGeneScore\n";
    foreach my $gene (sort {$icagesPrint{$b}{"score"} <=> $icagesPrint{$a}{"score"}} keys %icagesPrint){
        print GENES "$icagesPrint{$gene}{\"content\"}\n";
    }
}



