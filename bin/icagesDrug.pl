#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my ($rawInputFile, $icagesLocation, $prefix);
my (%biosystem, %neighbors, %activity, %onc, %sup, %icagesGenes);
# fda and clinical
# fda target
my (%fda, %clin, %fda_target, $fda_target_ref, $fdaRef, $clinRef);
my ($biosystemRef, $activityRef, $oncRef, $supRef, $icagesGenesRef, $neighborsRef, $supDrugRef, $oncDrugRef, $otherDrugRef);


######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

$rawInputFile = $ARGV[0];
$icagesLocation = $ARGV[1];
$prefix = $ARGV[2];
# add DGIdb database genes to reduce number of genes to query DGIdb
($biosystemRef, $activityRef, $oncRef, $supRef, $fdaRef, $clinRef, $fda_target_ref, $supDrugRef, $oncDrugRef, $otherDrugRef) = &loadDatabase($icagesLocation);

%biosystem = %{$biosystemRef};
%activity = %{$activityRef};
%onc = %{$oncRef};
%sup = %{$supRef};
%fda = %{$fdaRef};
%clin = %{$clinRef};
%fda_target = %{$fda_target_ref};
# %dgiGenes = %{$dgiGenesRef};
$icagesGenesRef = &getiCAGES($rawInputFile, $prefix);
%icagesGenes = %{$icagesGenesRef};

$neighborsRef = &getNeighbors(\%icagesGenes, \%biosystem, \%onc, \%sup);
%neighbors = %{$neighborsRef};
&getDrugs ($rawInputFile, $icagesLocation, \%neighbors, \%onc, \%sup, $prefix, \%fda_target, $supDrugRef, $oncDrugRef, $otherDrugRef);
&processDrugs($rawInputFile, \%neighbors, \%activity, $prefix , \%fda, \%clin); # add fda and clin

######################################################################################################################################
############################################################# subroutines ############################################################
######################################################################################################################################

sub loadDatabase {
    print "NOTICE: start loading Databases\n";
    my (%biosystem, %activity, %onc, %sup, %fda, %clin, %fda_target_ref, %supDrugDB, %oncDrugDB, %otherDrugDB);
    my ($icagesLocation, $DBLocation, $biosystemDB, $activityDB, $oncDB, $supDB, $fdaDB, $clinicalDB, $supDrugDB, $oncDrugDB, $otherDrugDB);
    $icagesLocation = shift;
    $DBLocation = $icagesLocation . "db/";
    $biosystemDB = $DBLocation . "biosystem.score";
    $activityDB = $DBLocation . "drug.score";
    $oncDB = $DBLocation . "oncogene.gene";
    $supDB = $DBLocation . "suppressor.gene";
    # fda and clinical trial, note that only one clinical trial would be provided, for more information, please visit myclinicaltrial.com
    $fdaDB = $DBLocation . "FDA_cancer.txt";
    $clinicalDB = $DBLocation . "ClinicalTrial.txt";
#    $dgiGenesDB = $DBLocation . "DGIdb.genes";
    $supDrugDB = $DBLocation . "suppressor.drug";
    $oncDrugDB = $DBLocation . "oncogene.drug";
    $otherDrugDB = $DBLocation . "othergene.drug";
    open(BIO, "$biosystemDB") or die "ERROR: cannot open $biosystemDB\n";
    open(ACT, "$activityDB") or die "ERROR: cannot open $activityDB\n";
    open(ONC, "$oncDB") or die "ERROR: cannot open $oncDB\n";
    open(SUP, "$supDB") or die "ERROR: cannot open $supDB\n";
    # fda and clinical trial
    open(FDA, "$fdaDB") or die "ERROR: cannot open $fdaDB\n";
    open(CLIN, "$clinicalDB") or die "ERROR: cannot open $clinicalDB\n";
 #   open(DGI, "$dgiGenesDB") or die "ERROR: cannot open $dgiGenesDB\n";
    open(SUPDRUG, "$supDrugDB") or die "ERROR: cannot open $supDrugDB\n";
    open(ONCDRUG, "$oncDrugDB") or die "ERROR: cannot open $oncDrugDB\n";
    open(OTHERDRUG, "$otherDrugDB") or die "ERROR: cannot open $otherDrugDB\n";

    while (<SUPDRUG>) {
        chomp;
        my $line = $_;
        my @line = split("\t", $_);
	next unless defined $line[0] and defined $line[1];
        $supDrugDB{$line[0]}{$line[1]} = $line;
#	print "$line[0]\t$line[1]\t$line\n";
    }

    while (<ONCDRUG>) {
        chomp;
        my $line = $_;
        my @line = split("\t", $_);
	next unless defined $line[0] and defined $line[1];
        $oncDrugDB{$line[0]}{$line[1]} = $line;
#	print "$line[0]\t$line[1]\t$line\n";
    }
    
    while (<OTHERDRUG>) {
        chomp;
        my $line = $_;
        my @line = split("\t", $_);
	next unless defined $line[0] and defined $line[1];
        $otherDrugDB{$line[0]}{$line[1]} = $line;
#	print "$line[0]\t$line[1]\t$line\n";
    }
    
  #  while(<DGI>){
   #     chomp;
    #    $dgiGenes{$_} = 1;
    #}
    
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
    
    while(<FDA>){
        chomp;
        my @line = split("\t", $_);
        # drugname: subtype,tradename
        $line[0] = "NA" unless defined $line[0];
        $line[1] = "NA" unless defined $line[1];
        $line[3] = "NA" unless defined $line[3];
        my $content = $line[1] . "," . $line[3];
        $fda{$line[0]} = $content;
        my @target = split(";", $line[2]);
        for(0..$#target){
            $fda_target{$target[$_]}{$line[0]} = 1;
        }
    }
    
    while(<CLIN>){
        chomp;
        my @line = split("\t", $_);
        # drugname: trialname,organization,phase,url
        $line[0] = "NA" unless defined $line[0];
        $line[1] = "NA"unless defined $line[1];
        $line[2] = "NA"unless defined $line[2];
        $line[3] = "NA" unless defined $line[3];
        $line[4] = "NA" unless defined $line[4];
        my $content = $line[1] . "," . $line[2] . "," . $line[3] . "," . $line[4];
        $clin{$line[0]} = $content;
    }
    
    close FDA;
    close CLIN;
    close SUP;
    close ONC;
    close ACT;
    close BIO;
    close SUPDRUG;
    close ONCDRUG;
    close OTHERDRUG;
    return (\%biosystem, \%activity, \%onc, \%sup, \%fda, \%clin, \%fda_target, \%supDrugDB, \%oncDrugDB, \%otherDrugDB);
}

sub getiCAGES{
    print "NOTICE: start process gene files from iCAGES layer two\n";
    my ($rawInputFile, $icagesGenes, $prefix);
    my %icagesGenes;
    $rawInputFile = shift;
    $prefix = shift;
    $icagesGenes = $rawInputFile . $prefix . ".annovar.icagesGenes.csv";
    open(GENES, "$icagesGenes") or die "ERROR: cannot open $icagesGenes\n";
    my $header = <GENES>;
    while(<GENES>){
        chomp;
        my @line = split(",", $_);
        $icagesGenes{$line[0]} = $line[5];
    }
    return \%icagesGenes;
}

sub getNeighbors{
    print "NOTICE: start getting top five neighbors for mutated genes\n";
    my (%icagesGenes, %biosystem, %neighbors, %onc, %sup);
    my ($icagesGenesRef, $biosystemRef, $oncRef, $supRef);
    my $index;
    $icagesGenesRef = shift;
    $biosystemRef = shift;
    $oncRef = shift;
    $supRef = shift;
#    $dgiGenesRef = shift;
    %icagesGenes = %{$icagesGenesRef};
    %biosystem = %{$biosystemRef};
    %onc = %{$oncRef};
    %sup = %{$supRef};

    foreach my $gene (sort keys %icagesGenes){
        $index = 0;
        $neighbors{$gene}{$gene}{"biosystem"} = 1;
        $neighbors{$gene}{$gene}{"icages"} = $icagesGenes{$gene};
        $neighbors{$gene}{$gene}{"product"} = $icagesGenes{$gene};
        foreach my $neighbor (sort { $biosystem{$b} <=> $biosystem{$a} }  keys %{$biosystem{$gene}}){
#	    if(exists $onc{$gene} or exists $sup{$gene}){
#		last if $index == 10;
#	    }else{
            last if $index == 5;
#	    }
            $index ++;
            $neighbors{$neighbor}{$gene}{"biosystem"} = $biosystem{$gene}{$neighbor};
            $neighbors{$neighbor}{$gene}{"icages"} = $icagesGenes{$gene};
            $neighbors{$neighbor}{$gene}{"product"} =  $icagesGenes{$gene} * $biosystem{$gene}{$neighbor};
        }
    }
    return \%neighbors;
}

sub getDrugs{
    print "NOTICE: start getting drugs for seed genes\n";
    my (%neighbors, %onc, %sup);
    my ($neighborsRef, $oncRef, $supRef);
    my (@seeds, @onc, @sup, @other);
    my ($onc, $sup, $other);
    $onc="";
    $sup="";
    $other="";
    # fda
    my $fda_target_ref;
    my %fda_target;

    my ($rawInputFile, $supFile, $oncFile, $otherFile, $icagesLocation, $callDgidb, $prefix);
    $rawInputFile = shift;
    $icagesLocation = shift;
    $neighborsRef = shift;
    $oncRef = shift;
    $supRef = shift;
    $prefix = shift;
    $fda_target_ref = shift;
    # three kinds of drugs
    my $supDrugRef = shift;
    my $oncDrugRef = shift;
    my $otherDrugRef = shift;
    my %supDrug;
    my %oncDrug;
    my %otherDrug;
    %supDrug = %{$supDrugRef};
    %oncDrug = %{$oncDrugRef};
    %otherDrug = %{$otherDrugRef};
    
    %neighbors = %{$neighborsRef};
    %onc = %{$oncRef};
    %sup = %{$supRef};
    %fda_target = %{$fda_target_ref};
    @seeds = keys %neighbors;
    $callDgidb = $icagesLocation . "bin/DGIdb/getDrugList.pl";
    $supFile = $rawInputFile . $prefix . ".suppressor.drug";
    $oncFile = $rawInputFile . $prefix.".oncogene.drug";
    $otherFile = $rawInputFile . $prefix. ".other.drug";

    # find FDA drugs for genes in case DGIdb missed it
    my %all_genes;
    for(0..$#seeds){
#        print "$seeds[$_]\n";
        if(exists $sup{$seeds[$_]} and $seeds[$_] =~ /[a-zA-Z0-9]+/){
#            if(exists $dgiGenes{$seeds[$_]}){
                push @sup, $seeds[$_];
#            }
            $all_genes{$seeds[$_]} = 1;
        }elsif(exists $onc{$seeds[$_]} and $seeds[$_] =~ /[a-zA-Z0-9]+/){
#            if(exists $dgiGenes{$seeds[$_]}){
                push @onc, $seeds[$_];
#            }
            $all_genes{$seeds[$_]} = 1;
        }else{
            if($seeds[$_] =~ /[a-zA-Z0-9]+/){
#                if(exists $dgiGenes{$seeds[$_]}){
                    push @other, $seeds[$_];
 #               }
            $all_genes{$seeds[$_]} = 1;
            }
        }
    }
    $sup = join(",", @sup);
    $onc = join(",", @onc);
    $other = join(",", @other);
   
    open(SUP, ">$supFile") or die;
    open(ONC, ">$oncFile") or die;
    open(OTHER, ">$otherFile") or die;
    
    
    if($sup ne "" ){
#	print "iCAGES: $callDgidb --genes='$sup' --interaction_type='activator,other/unknown,n/a,inducer,positive allosteric modulator,potentiator,stimulator' --source_trust_levels='Expert curated' --output='$supFile'";
        for(0..$#sup) {
#	    print "$sup[$_]\n";
           if (exists $supDrug{$sup[$_]}) {
                foreach my $key (sort keys %{$supDrug{$sup[$_]}}) {
                    print SUP "$supDrug{$sup[$_]}{$key}\n";
                }
            }
        }
        
        # !system("$callDgidb --genes='$sup' --interaction_type='activator,other/unknown,n/a,inducer,positive allosteric modulator,potentiator,stimulator' --source_trust_levels='Expert curated' --output='$supFile'") or warn "ERROR: cannot gt drugs\n$callDgidb --genes='$sup' --interaction_type='activator,other/unknown,n/a,inducer,positive allosteric modulator,potentiator,stimulator' --source_trust_levels='Expert curated' --output='$supFile'";
    }
    if($onc ne ""){
#	print "iCAGES: $callDgidb --genes='$onc' --interaction_type='agonist,antisense,competitive,immunotherapy,inhibitory allosteric modulator,inverse agonist,negative modulator,partial agonist,partial antagonist,vaccine,inhibitor,suppressor,antibody,antagonist,blocker,other/unknown,n/a' --source_trust_levels='Expert curated' --output='$oncFile'";
        #!system("$callDgidb --genes='$onc' --interaction_type='agonist,antisense,competitive,immunotherapy,inhibitory allosteric modulator,inverse agonist,negative modulator,partial agonist,partial antagonist,vaccine,inhibitor,suppressor,antibody,antagonist,blocker,other/unknown,n/a' --source_trust_levels='Expert curated' --output='$oncFile'") or warn "ERROR: cannot get drugs\n$callDgidb --genes='$onc' --interaction_type='agonist,antisense,competitive,immunotherapy,inhibitory allosteric modulator,inverse agonist,negative modulator,partial agonist,partial antagonist,vaccine,inhibitor,suppressor,antibody,antagonist,blocker,other/unknown,n/a' --source_trust_levels='Expert curated' --output='$oncFile'\n";
        
        for(0..$#onc) {
#	    print "$onc[$_]\n";
            if (exists $oncDrug{$onc[$_]}) {
                foreach my $key (sort keys %{$oncDrug{$onc[$_]}}) {
                    print ONC "$oncDrug{$onc[$_]}{$key}\n";
                }
            }
        }
    }
    if($other ne ""){
        #	print "$callDgidb --genes='$other' --source_trust_levels='Expert curated' --output='$otherFile'";
        
        #!system("$callDgidb --genes='$other' --source_trust_levels='Expert curated' --output='$otherFile'") or warn "ERROR: cannot get drugs\n$callDgidb --genes='$other' --source_trust_levels='Expert curated' --output='$otherFile'\n";
        
        for(0..$#other) {
#	    print "$other[$_]\n";
            if (exists $otherDrug{$other[$_]}) {
                foreach my $key (sort keys %{$otherDrug{$other[$_]}}) {
                    print OTHER "$otherDrug{$other[$_]}{$key}\n";
                }
            }
        }
        
    }
    
    # create an output file
    my $fdaDrug = $rawInputFile . $prefix. ".fda.drug";
    open(FDADRUG, ">$fdaDrug") or die "ERROR: cannot create $fdaDrug for outputing FDA drugs\n";
    foreach my $gene (sort keys %all_genes){
	if(exists $fda_target{$gene}){
	    foreach my $drug (sort keys %{$fda_target{$gene}}){
		print FDADRUG "$gene\t$drug\n";
	    }
	}
    }
}


sub processDrugs{
    print "NOTICE: start processing drugs from DGIdb\n";
    my ($rawInputFile, $matchFile, $allDrugs, $icagesDrugs, $prefix);
    my (%neighbors, %activity, %icagesDrug, %icagesPrint);
    my ($neighborsRef, $activityRef);
    my ($oncDrugFile, $supDrugFile, $otherDrugFile);
    # fda and clin
    my ($fdaRef, $clinRef);
    my (%fda, %clin);
    # fda drug list
    my $FDADrugFile;

    $rawInputFile = shift;
    $neighborsRef = shift;
    $activityRef = shift;
    $prefix = shift;
    $fdaRef = shift;
    $clinRef = shift;
    %fda = %{$fdaRef};
    %clin = %{$clinRef};
    %neighbors = %{$neighborsRef};
    %activity = %{$activityRef};
    $matchFile = $rawInputFile . $prefix . ".*.drug";
    $oncDrugFile = $rawInputFile . $prefix . ".oncogene.drug";
    $supDrugFile = $rawInputFile . $prefix . ".suppressor.drug";
    $otherDrugFile = $rawInputFile . $prefix. ".other.drug";
    $FDADrugFile = $rawInputFile . $prefix. ".fda.drug";
    $allDrugs = $rawInputFile . $prefix . ".drug.all";
    $icagesDrugs = $rawInputFile . $prefix . ".annovar.icagesDrugs.csv";
    if(! -e $oncDrugFile){
	!system("touch $oncDrugFile") or die "ERROR: cannot create $oncDrugFile\n";
    }
    if(! -e $supDrugFile){
	!system("touch $supDrugFile") or die "ERROR: cannot create $supDrugFile\n";
    }
    if(! -e $otherDrugFile){
	!system("touch $otherDrugFile") or die "ERROR: cannot create $otherDrugFile\n";
    }
    if(! -e $FDADrugFile){
        !system("touch $FDADrugFile") or die "ERROR: cannot create $FDADrugFile\n";
    }
    !system("cat $oncDrugFile $supDrugFile $otherDrugFile $FDADrugFile > $allDrugs") or die "ERROR: cannot create an empty drug file\n";
    open(DRUG, "$allDrugs") or die "ERROR: cannot open drug file $allDrugs\n";
    open(OUT, ">$icagesDrugs") or die "ERROR: cannot open $icagesDrugs\n";
    while(<DRUG>){
        chomp;
        my @line = split("\t", $_);
        next unless defined $line[1];
        my $neighbor = $line[0];
        my $index = 0;
        foreach my $target (sort { $neighbors{$neighbor}{$b}{"product"} <=> $neighbors{$neighbor}{$a}{"product"} } keys %{$neighbors{$neighbor}}){
            last if $index == 1;
            if(exists $icagesDrug{$line[1]}{$neighbor}){
                if($neighbors{$neighbor}{$target}{"product"} > $icagesDrug{$line[1]}{$neighbor}{$target}{"biosystem"} * $icagesDrug{$line[1]}{$neighbor}{$target}{"icages"}){
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"biosystem"} = $neighbors{$neighbor}{$target}{"biosystem"};
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"icages"} = $neighbors{$neighbor}{$target}{"icages"} ;
		    if(exists $activity{$line[1]}){
                        $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = $activity{$line[1]};
                    }else{
                        $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = 0;
                    }
#		    $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = 1 if exists $fda{$line[1]} ;
#		    $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = max(0.5, $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"}) if  exists $clin{$line[1]};
                }
            }else{
                $icagesDrug{$line[1]}{$neighbor}{$target}{"biosystem"} = $neighbors{$neighbor}{$target}{"biosystem"};
                $icagesDrug{$line[1]}{$neighbor}{$target}{"icages"} = $neighbors{$neighbor}{$target}{"icages"} ;
                if(exists $activity{$line[1]}){
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = $activity{$line[1]};
                }else{
                    $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = 0;
                }
#		$icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = 1 if exists $fda{$line[1]} ;
#		$icagesDrug{$line[1]}{$neighbor}{$target}{"activity"} = max(0.5, $icagesDrug{$line[1]}{$neighbor}{$target}{"activity"}) if exists $clin{$line[1]};
		
            }
            $index ++;
        }
    }
    
    ##### count drug
    my $drugCount = 0;
    my $gooddrugCount = 0;
    foreach my $drug (sort keys %icagesDrug){
        foreach my $neighbor (sort keys %{$icagesDrug{$drug}}){
            foreach my $final (sort keys %{$icagesDrug{$drug}{$neighbor}}){
                my $icagesDrug = $icagesDrug{$drug}{$neighbor}{$final}{"biosystem"} * $icagesDrug{$drug}{$neighbor}{$final}{"icages"} * $icagesDrug{$drug}{$neighbor}{$final}{"activity"};
		my $tier = 3;
		if(exists $fda{$drug}){
		    $tier = 1;
		}elsif(exists $clin{$drug}){
		    $tier = 2;
		}
                $icagesPrint{$tier}{$drug}{"score"} = $icagesDrug;

		if( $icagesDrug{$drug}{$neighbor}{$final}{"activity"} == 0){
		    $icagesDrug{$drug}{$neighbor}{$final}{"activity"} = "NA";
		}
		if($icagesDrug == 0){
		    $icagesDrug{$drug}{$neighbor}{$final}{"drug"} ="NA";
		}
                $icagesPrint{$tier}{$drug}{"content"} = "$drug,$final,$neighbor,$icagesDrug{$drug}{$neighbor}{$final}{\"icages\"},$icagesDrug{$drug}{$neighbor}{$final}{\"biosystem\"},$icagesDrug{$drug}{$neighbor}{$final}{\"activity\"},$icagesDrug,$tier";
            }
        }
    }
    print OUT "drugName,finalTarget,directTarget,iCAGESGeneScore,maxBioSystemsScore,maxActivityScore,icagesDrugScore,tier,FDA_approvedSubtype,FDA_activeIngredient,CLT_name,CLT_organization,CLT_phase,CLT_url\n";
    foreach my $tier (sort {$a <=> $b} keys %icagesPrint){
    foreach my $drug (sort {$icagesPrint{$tier}{$b}{"score"} <=> $icagesPrint{$tier}{$a}{"score"}} keys %{$icagesPrint{$tier}}){
        $drugCount ++;
        $gooddrugCount ++ if $icagesPrint{$tier}{$drug}{"score"} > 0.5;
	# check fda and clinical trial
	my $printContent = $icagesPrint{$tier}{$drug}{"content"};
	if(exists $fda{$drug}){
	    $printContent .= "," . $fda{$drug};
	}else{
	    $printContent .= ",NA,NA";
	}
	if(exists $clin{$drug}){
	    $printContent .= "," . $clin{$drug};
	}else{
	    $printContent .= ",NA,NA,NA,NA";
	}
        print OUT "$printContent\n";
    }
    }
    close OUT;
    close DRUG;
    
    my $logFile = $rawInputFile . $prefix . ".annovar.icages.log";
    open(LOG, ">>$logFile") or die "iCAGES: cannot open file $logFile\n";
    print LOG "########### iCAGES Drug Summary ###########\n";
    print LOG "## basic information\n";
    print LOG "Total: $drugCount\n";
    print LOG "Good drug (iCAGES drug score >= 0.5): $gooddrugCount\n";
}







