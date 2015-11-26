#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);
use Pod::Usage;
use Getopt::Long;
use JSON;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my $rawInputFile = $ARGV[0];
my $icagesLocation = $ARGV[1];
my $prefix = $ARGV[2];
my ($icagesMutationsRef, $icagesGenesRef, $icagesDrugsRef, $logInformationRef);
my $json;
my $nowString;
$nowString = localtime();

######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

($icagesMutationsRef, $icagesGenesRef, $icagesDrugsRef, $logInformationRef) = &loadResult($rawInputFile, $prefix);
$json = &createJson($icagesMutationsRef, $icagesGenesRef, $icagesDrugsRef, $logInformationRef);
&printJson($rawInputFile, $json, $prefix);

######################################################################################################################################
############################################################# subroutines ############################################################
######################################################################################################################################

sub loadResult {
    print "NOTICE: start loading three output files from iCAGES\n";
    my ($rawInputFile, $prefix, $icagesMutationsFile, $icagesGenesFile, $icagesDrugsFile);
    my ($icagesMutationsRef, $icagesGenesRef, $icagesDrugsRef);
    my (%icagesMutations, %icagesGenes, %icagesDrugs);
    my ($missenseCount, $noncodingCount, $structuralVariationCount);
    my ($geneCount, $driverCount, $cgcCount, $keggCount, $drugCount);
    my %logInformation;
    $rawInputFile = shift;
    $prefix = shift;
    $icagesMutationsFile = $rawInputFile . $prefix. ".annovar.icagesMutations.csv";
    $icagesGenesFile = $rawInputFile . $prefix. ".annovar.icagesGenes.csv";
    $icagesDrugsFile = $rawInputFile . $prefix. ".annovar.icagesDrugs.csv";
    ($icagesMutationsRef, $missenseCount, $noncodingCount, $structuralVariationCount) = &loadMutations($icagesMutationsFile);
    ($icagesGenesRef, $geneCount, $driverCount, $cgcCount, $keggCount) = &loadGenes($icagesGenesFile);
    ($icagesDrugsRef, $drugCount) = &loadDrugs($icagesDrugsFile);
    %icagesMutations = %{$icagesMutationsRef};
    %icagesGenes = %{$icagesGenesRef};
    %icagesDrugs = %{$icagesDrugsRef};
    %logInformation = ("Gene_count" => $geneCount, "Driver_count" => $driverCount, "CGC_count" => $cgcCount, "KEGG_count" => $keggCount, "Missense_count" => $missenseCount, "Noncoding_count" => $noncodingCount, "Structural_variation_count" => $structuralVariationCount, "Drug_count" => $drugCount);
    return (\%icagesMutations, \%icagesGenes, \%icagesDrugs, \%logInformation);
}

sub createJson {
    print "NOTICE: start creating JSON file\n";
    my ($icagesMutationsRef, $icagesGenesRef, $icagesDrugsRef, $logInformationRef);
    my (%icagesMutations, %icagesGenes, %icagesDrugs);
    my %json;
    my $json;
    $icagesMutationsRef = shift;
    $icagesGenesRef = shift;
    $icagesDrugsRef = shift;
    $logInformationRef = shift;
    %icagesMutations = %{$icagesMutationsRef};
    %icagesGenes = %{$icagesGenesRef};
    %icagesDrugs = %{$icagesDrugsRef};
    foreach my $gene (sort keys %icagesGenes){
        $icagesGenes{$gene}{"Mutation"} = $icagesMutations{$gene};
        $icagesGenes{$gene}{"Children"} = $icagesDrugs{$gene};
        push @{$json{"Output"}}, $icagesGenes{$gene};
    }
    $json{"Log"} = $logInformationRef;
    $json = encode_json \%json;
    return $json;
}

sub printJson {
    my ($rawInputFile, $json, $prefix, $icagesJsonFile);
    $rawInputFile = shift;
    $json = shift;
    $prefix = shift;
    $icagesJsonFile = $rawInputFile . $prefix . ".icages.json";
    open(OUT, ">$icagesJsonFile") or die ;
    print OUT "$json\n";
    close OUT;
    print "NOTICE: end runing iCAGES packge at $nowString\n";
}

sub loadMutations {
    print "NOTICE: start processing icagesMutations.csv file to generate JSON file\n";
    my ($icagesMutationsFile, $title);
    my ($missenseCount, $noncodingCount, $structuralVariationCount);
    my %icagesMutations;
    $icagesMutationsFile = shift;
    open(MUT, "$icagesMutationsFile") or die;
    $title = <MUT>;
    $missenseCount = 0;
    $noncodingCount = 0;
    $structuralVariationCount = 0;
    while(<MUT>){
        chomp;
        my @line = split(",", $_);
        my $geneName = $line[0];
        my $chrmosomeNumber = $line[1];
        my $start = $line[2];
        my $end = $line[3];
        my $reference = $line[4];
        my $alternative = $line[5];
        my $category = $line[6];
        if($category eq "point coding"){
            $missenseCount ++;
        }elsif($category eq  "point noncoding"){
            $noncodingCount ++;
        }elsif($category eq  "structural variation"){
            $structuralVariationCount ++;
        }
        my $mutationSyntax = $line[7];
        my $proteinSyntax = $line[8];
        my $scoreCategory = $line[9];
        my $mutationScore = $line[10] eq "NA" ? "NA" : sprintf("%.2f", $line[10]);
	
        my %anonymousHash = ("Chromosome" => $chrmosomeNumber, "Start_position" => $start, "End_position" => $end, "Reference_allele" => $reference, "Alternative_allele" => $alternative, "Mutation_syntax" => $mutationSyntax, "Protein_syntax" => $proteinSyntax, "Mutation_category" => $category, "Score_category" => $scoreCategory, "Driver_mutation_score" => $mutationScore);
        push @{$icagesMutations{$geneName}}, \%anonymousHash;
    }
    close MUT;
    return (\%icagesMutations, $missenseCount, $noncodingCount, $structuralVariationCount);
}

sub loadGenes {
    print "NOTICE: start processing icagesGenes.csv file to generate JSON file\n";
    my ($icagesGenesFile, $title);
    my ($geneCount, $driverCount, $cgcCount, $keggCount);
    $geneCount = 0;
    $driverCount = 0;
    $cgcCount = 0;
    $keggCount = 0;
    $icagesGenesFile = shift;
    my $keggRef;
    my %kegg;
    my %icagesGenes;
    my $DBLocation = $icagesLocation . "db/";
    $keggRef = loadDatabase($DBLocation);
    %kegg = %{$keggRef};
    open(GENE, "$icagesGenesFile") or die;
    $title = <GENE>;
    # limit the number of maximum genes of 20 
    while(<GENE>){
        chomp;
	last if $geneCount == 49;
        $geneCount ++;
        my @line = split(",", $_);
        my $geneName = $line[0];
        my $category = $line[1];
	my $phenolyzer = $line[5] eq "NA" ? "NA" : sprintf("%.2f", $line[5]);
	my $icagesGeneScore = $line[6] eq "NA" ? "NA" : sprintf("%.2f", $line[6]);
        my $url;
        if ($category eq "Cancer Gene Census"){
            $cgcCount ++;
            $url = "http://cancer.sanger.ac.uk/cosmic/gene/overview?ln=" . $geneName;
        }elsif($category eq "KEGG Cancer Pathway"){
            $keggCount ++;
            $url = "http://www.genome.jp/dbget-bin/www_bget?hsa:" . $kegg{$geneName};
        }else{
            $url = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=". $geneName. "&search=96a88d95c4a24ffc2ac1129a92af7b02";
        }
        my %anonymousHash =  ("Gene_url" => $url, "Name" => $geneName, "Category" => $category, "Phenolyzer_score" => $phenolyzer, "iCAGES_gene_score" => $icagesGeneScore);
        $icagesGenes{$geneName} = \%anonymousHash;
    }
    my @sortediCAGESGeneScore = sort {$icagesGenes{$b}{"iCAGES_gene_score"} <=> $icagesGenes{$a}{"iCAGES_gene_score"}} keys %icagesGenes;
    my $percentile = int(($#sortediCAGESGeneScore + 1) * 0.20);
    my $criticalValue = $icagesGenes{$sortediCAGESGeneScore[$percentile-1]}{"iCAGES_gene_score"};
    foreach my $gene (sort keys %icagesGenes){
        if($icagesGenes{$gene}{"iCAGES_gene_score"} >= $criticalValue){
            $icagesGenes{$gene}{"Driver"} = "TRUE";
            $driverCount ++;
        }else{
            $icagesGenes{$gene}{"Driver"} = "FALSE";
        }
    }
    close GENE;
    return (\%icagesGenes, $geneCount, $driverCount, $cgcCount, $keggCount);
}



sub loadDrugs {
    print "NOTICE: start processing icagesGenes.csv file to generate JSON file\n";
    my ($icagesDrugsFile, $title);
    my %icagesDrugs;
    my $drug_count;
    $drug_count = 0;
    $icagesDrugsFile = shift;
    open(DRUG, "$icagesDrugsFile") or die;
    $title = <DRUG>;
    while(<DRUG>){
        $drug_count ++;
        chomp;
        my @line = split(",", $_);
        my $drugName = $line[0];
        my $finalTarget = $line[1];
        my $directTarget = $line[2];
	my $maxBioSystemsScore = $line[4] eq "NA" ? "NA" : sprintf("%.2f", $line[4]);
	my $maxActivityScore = $line[5] eq "NA" ? "NA" : sprintf("%.2f", $line[5]);
	my $icagesDrugScore = $line[6] eq "NA" ? "NA" : sprintf("%.2f", $line[6]);
	# add fda(7,8) and clinical trial(9,10,11,12)
	my $FDA_tag = ( $line[7] eq "NA" and $line[8] eq "NA" ) ? "FALSE" : "TRUE" ;
	my $CT_tag = ( $line[9] eq "NA" and $line[10] eq "NA" and $line[11] eq "NA" and $line[12] eq "NA" ) ? "FALSE" : "TRUE" ;
	my (%fda, %ct);
        my %anonymousHash;
	if($FDA_tag eq "FALSE" and $CT_tag eq "FALSE"){
	    %anonymousHash = ("Drug_name" => $drugName, "Final_target_gene" => $finalTarget, "Direct_target_gene" => $directTarget, "BioSystems_probability" => $maxBioSystemsScore, "PubChem_active_probability" => $maxActivityScore, "iCAGES_drug_score" => $icagesDrugScore, "Target_mutation_tag" => "FALSE" , "FDA_tag" => "FALSE", "CT_tag" => "FALSE");
	}elsif($FDA_tag eq "TRUE" and $CT_tag eq "FALSE"){
	    %fda = ("Status" => $line[7], "Active_ingredient" => $line[8]);
	    %anonymousHash = ("Drug_name" => $drugName, "Final_target_gene" => $finalTarget, "Direct_target_gene" => $directTarget, "BioSystems_probability" => $maxBioSystemsScore, "PubChem_active_probability" => $maxActivityScore, "iCAGES_drug_score" => $icagesDrugScore, "Target_mutation_tag" => "FALSE" , "FDA_tag" => "TRUE", "CT_tag" => "FALSE", "FDA_Info" => \%fda);
	}elsif($CT_tag eq "TRUE" and $FDA_tag eq "FALSE"){
	    %ct = ("Name" => $line[9], "Organization" => $line[10], "Phase" => $line[11] , "URL" => $line[12]);
	    %anonymousHash = ("Drug_name" => $drugName, "Final_target_gene" => $finalTarget, "Direct_target_gene" => $directTarget, "BioSystems_probability" => $maxBioSystemsScore, "PubChem_active_probability" => $maxActivityScore, "iCAGES_drug_score" => $icagesDrugScore, "Target_mutation_tag" => "FALSE" , "FDA_tag" => "FALSE", "CT_tag" => "TRUE", "CT_Children"=> [\%ct]);
	}else{
	    %fda = ("Status" => $line[7], "Active_ingredient" => $line[8]);
	    %ct= ("Name" => $line[9], "Organization" => $line[10], "Phase" => $line[11] , "URL" => $line[12]);
	    %anonymousHash = ("Drug_name" => $drugName, "Final_target_gene" => $finalTarget, "Direct_target_gene" => $directTarget, "BioSystems_probability" => $maxBioSystemsScore, "PubChem_active_probability" => $maxActivityScore, "iCAGES_drug_score" => $icagesDrugScore, "Target_mutation_tag" => "FALSE" , "FDA_tag" => "TRUE", "CT_tag" => "TRUE", "FDA_Info" => \%fda, "CT_Children" => [\%ct]);
	}
	
        push @{$icagesDrugs{$finalTarget}}, \%anonymousHash;
    }
    close DRUG;
    return (\%icagesDrugs, $drug_count);
}


sub loadDatabase {
    print "NOTICE: start loading CGC, KEGG cancer pathway databases";
    my ($DBLocation, $cgcFile, $keggFile);
    my (%cgc, %kegg);
    $DBLocation = shift;
    $keggFile = $DBLocation . "kegg.gene";
    open(KEGG, "$keggFile") or die ;
    while(<KEGG>){
        chomp;
        my @line = split;
        $kegg{$line[0]} = $line[1];
    }
    close KEGG;
    return (\%kegg);
}


