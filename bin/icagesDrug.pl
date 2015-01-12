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
my (%biosystem, %activity, %onc, %sup, %icagesGenes);

######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

$rawInputFile = $ARGV[0];
$icagesLocation = $ARGV[1];
(%biosystem, %activity, %onc, %sup) = &loadDatabase($icagesLocation);
%icagesGenes = &getiCAGES($rawInputFile);
%neighbors = &getNeighbors(\%icagesGenes, \%biosystem);
&getDrugs ($rawInputFile, $icagesLocation, \%neighbors, \%onc, \%sup);
&processDrugs($rawInputFile, \%neighbors, \%activity);

######################################################################################################################################
############################################################# subroutines ############################################################
######################################################################################################################################

sub loadDatabase {
    print "NOTICE: start loading Databases\n";
    my (%biosystem, %activity, %onc, %sup);
    my ($icagesLocation, $DBLocation, $biosystemDB, $activityDB, $oncDB, $supDB);
    $icagesLocation = shift;
    $DBLocation = $icagesLocation . "db/";
    $biosystemDB = $DBLocation . "biosystem.score";
    $activityDB = $DBLocation . "drug.score";
    $oncDB = $DBLocation . "oncogene.gene";
    $supDB = $DBLocation . "suppressor.gene";
    open(BIO, "$biosystemDB") or die "ERROR: cannot open $biosystemDB\n";
    open(ACT, "$activityDB") or die "ERROR: cannot open $activityDB\n";
    open(ONC, "$oncDB") or die "ERROR: cannot open $oncDB\n";
    open(SUP, "$supDB") or die "ERROR: cannot open $supDB\n";
    while(<BIO>){
        chomp;
        my @line = split("\t", $_);
        $biosystem{$line[0]}{$line[1]} = $line[2];
    }
    while(<ACT>){
        chomp;
        my @line = split("\t", $_);
        $activity{$line[0]} = $line[1];
    }
    while(<ONC>){
        chomp;
        my @line = split("\t", $_);
        $onc{$line[0]} = 1;

    }
    while(<SUP>){
        chomp;
        my @line = split("\t", $_);
        $sup{$line[0]} = 1;
        
    }
    close SUP;
    close ONC;
    close ACT;
    close BIO;
    return (\%biosystem, \%activity, \%onc, \%sup);
}

sub getiCAGES{
    print "NOTICE: start process gene files from iCAGES layer two\n";
    my ($rawInputFile, $icagesGenes);
    my %icagesGenes;
    $rawInputFile = shift;
    $icagesGenes = $rawInputFile . ".annovar.icagesGenes.csv";
    open(GENES, "$icagesGenes") or die "ERROR: cannot open $icagesGenes\n";
    while(<GENES>){
        chomp;
        my @line = split(",", $_);
        $icagesGene{$line[0]} = $line[6];
    }
    return \%icagesGenes;
}

sub getNeighbors{
    print "NOTICE: start getting top five neighbors for mutated genes\n";
    my (%icagesGenes, %biosystem, %neighbors);
    my $index;
    %icagesGenes = shift;
    %biosystem = shift;
    foreach my $gene (sort keys %icagesGenes){
        $index = 0;
        $neighbors{$gene}{$gene}{"biosystem"} = 1;
        $neighbors{$gene}{$gene}{"icages"} = $icagesGenes{$gene};
        foreach my $neighbor (sort { $biosystem{$b} <=> $biosystem{$a} }  keys %{$biosystem{$gene}}){
            last if $index == 5;
            $index ++;
            $neighbors{$neighbor}{$gene}{"biosystem"} = $biosystem{$gene}{$neighbor};
            $neighbors{$neighbor}{$gene}{"icages"} = $icagesGenes{$gene};
        }
    }
    return \%neighbors;
}

sub getDrugs{
    print "NOTICE: start getting drugs for seed genes\n";
    my (%neighbors, %onc, %sup);
    my (@seeds, @onc, @sup, @other);
    my ($onc, $sup, $other);
    my ($rawInputFile, $supFile, $oncFile, $otherFile, $icagesLocation, $callDgidb);
    $rawInputFile = shift;
    $icagesLocation = shift;
    %neighbors = shift;
    %onc = shift;
    %sup = shift;
    @seeds = keys %neighbors;
    $callDgidb = $icagesLocation . "bin/DGIdb/getDrugList.pl";
    $supFile = $rawInputFile . ".suppressor.drug";
    $oncFile = $rawInputFile . ".oncogene.drug";
    $otherFile = $rawInputFile . ".other.drug";
    for(0..$#seeds){
        if(exists $sup{$seeds[$_]}){
            push @sup, $seeds[$_];
        }elsif(exists $onc{$seeds[$_]}){
            push @onc, $seeds[$_];
        }else{
            push @other, $seeds[$_];
        }
    }
    $sup = join(",", @sup);
    $onc = join(",", @onc);
    $other = join(",", @other);
    if($suppressors ne ""){
        !system("$callDgidb --genes='$sup' --interaction_type='activator,other/unknown,n/a,inducer,stimulator' --source_trust_levels='Expert curated' --output='$supFile'") or die "ERROR: cannot get drugs\n";
    }
    if($oncogenes ne ""){
        !system("$callDgidb --genes='$onc' --interaction_type='inhibitor,suppressor,antibody,antagonist,blocker,other/unknown,n/a' --source_trust_levels='Expert curated' --output='$oncFile'") or die "ERROR: cannot get drugs\n";
    }
    if($othergenes ne ""){
        !system("$callDgidb --genes='$other' --interaction_type='inhibitor,suppressor,antibody,antagonist,blocker,activator,other/unknown,n/a,inducer,stimulator' --source_trust_levels='Expert curated' --output='$otherFile'") or die "ERROR: cannot get drugs\n";
    }
}


sub processDrugs{
    print "NOTICE: start processing drugs from DGIdb\n";
    my ($rawInputFile, $matchFile, $allDrugs, $icagesDrugs);
    my %neighbors;
    my %activity;
    my %icagesDrug;
    $rawInputFile = shift;
    $matchFile = $rawInputFile . ".*.drug";
    $allDrugs = $rawInputFile . ".all.drug";
    $icagesDrugs = $rawInputFile . ".icagesDrugs.csv";
    !system("cat $matchFile > $allDrugs") or die "ERROR: cannot concatenate drug files\n";
    open(DRUG, "$allDrugs") or die "ERROR: cannot open drug file $allDrugs\n";
    open(OUT, "$icagesDrugs") or die "ERROR: cannot open $icagesDrugs\n";
    while(<DRUG>){
        chomp;
        my @line = split("\t", $_);
        my $neighbor = $line[0]
        my $index = 0;
        foreach my $target (sort { $neighbor{$neighbor}{$b}{"biosystem"}*$neighbor{$neighbor}{$b}{"icages"} <=> $neighbor{$neighbor}{$a}{"biosystem"}*$neighbor{$neighbor}{$a}{"icages"} }){
            last if $index == 1;
            if(exists $drug{$line[1]}{$neighbor}){
                if($neighbor{$neighbor}{$target}{"biosystem"}*$neighbor{$target}{"icages"} > $drug{$line[1]}{$neighbor}{$target}{"biosystem"} * $drug{$line[1]}{$neighbor}{$target}{"icages"}){
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"biosystem"} = $neighbor{$neighbor}{$target}{"biosystem"};
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"icages"} = $neighbor{$neighbor}{$target}{"icages"} ;
                    if(exists $activity{$line[1]}){
                        $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = $activity{$line[1]};
                    }else{
                        $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = 0;
                    }
                }
            }else{
                $icagesDrug{$line[1]}{$neighbor}{$target}{"biosystem"} = $neighbor{$neighbor}{$target}{"biosystem"};
                $icagesDrug{$line[1]}{$neighbor}{$target}{"icages"} = $neighbor{$neighbor}{$target}{"icages"} ;
                if(exists $activity{$line[1]}){
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = $activity{$line[1]};
                }else{
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = 0;
                }
            }
            $index ++;
        }
    }
    foreach my $drug (sort keys %icagesDrug){
        foreach my $neighbor (sort keys %{$icagesDrug{$drug}}){
            foreach my $final (sort keys %{$icagesDrug{$drug}{$neighbor}}){
                my $icagesDrug = $icagesDrug{$drug}{$neighbor}{"biosystem"} * $icagesDrug{$drug}{$neighbor}{"icages"} * $icagesDrug{$drug}{$neighbor}{"activity"};
                print OUT "$drug,$final,$neighbor,$icagesDrug{$drug}{$neighbor}{\"icages\"},$icagesDrug{$drug}{$neighbor}{\"biosystem\"},$icagesDrug{$drug}{$neighbor}{\"activity\"},$icagesDrug\n";
            }
        }
    }
    close DRUG;
}







