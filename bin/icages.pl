#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);	
use Pod::Usage;
use Getopt::Long;


######################################################################################################################################
########################################################### This is iCAGES! ##########################################################
######################################################################################################################################

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my ($inputFile, $outputFile, $snpFile, $cnvFile, $gene, $icages, $funseq, $cnv, $annovarLog);       #ANNOVAR files
my ($annovarCommand, $dgidbCommand, $callDgidb, $callConvertToAnnovar, $callAnnovar);               #ANNOVAR commands
my ($icagesDB, $icagesIndex, $funseqDB, $funseqIndex, $refGeneDB, $refGeneIndex, $cnvDB);           #ANNOVAR DB files: iCAGES score (index), refGene (fasta), dbSNP
my ($formatCheckFirstLine);                                                                         #ANNOVAR checking: format, missense mutations

my ($huiDB, $cgcDB, $keggDB, $suppressorDB, $oncogeneDB, $zscoreDB, $biosystemDB);                  #iCAGES DB files: Phenolyzer score, cgc genes, kegg cancer pathway genes
my ($help, $percent, $manual, $downloadDB, $id);                                                    #iCAGES arguments
my ($txt, $drug, $log);                                                                             #iCAGES output files
my @log;                                                                                            #iCAGES log content
my (@oncogenes, @suppressors, @othergenes, @relatedgenes);


my ($suppressordrugs, $oncogenedrugs, $otherdrugs, $relateddrugs);                                  #DGIdb temp files
my @drugs;

my $perlCommand;                                                                                    #SYSTEM commands
my ($icagesLocation, $funseqLocation, $inputLocation, $outputLocation, $DBLocation, $tempLocation, $logLocation);    #SYSTEM locations
my $nowString = localtime;                                                                          #SYSTEM local time

my %radialSVM;                                                                                      #HASH mutation->radialSVM
my %icages;                                                                                         #HASH geneName->[0]icages(total), [1]icages(highest), [2]cnv(highest)
my %hui;                                                                                            #HASH geneName->hui
my %kegg;                                                                                           #HASH geneName->kegg
my %cgc;                                                                                            #HASH geneName->1
my %icages_txt;                                                                                     #HASH geneName->mutation->information
my %cnv;                                                                                            #HASH mutation->1
my %sup;
my %onc;
my %zscore;
my %drug;
my %driver;                                                                                         #HASH gene->1
my %related;                                                                                        #HASH driver->related->score
my %funseq;

my ($minHui, $maxHui);                                                                              #MIN/MAX for Phenolyzer scores
my ($minGene, $maxGene);                                                                            #MIN/MAX for iCAGES score for all genes (for normalization purpose)
my ($minzscore, $maxzscore);

my $geneNumber;                                                                                     #JSON number of genes
my $iGene;
my $iDrug;

my (@radialSVMUnsort, @radialSVMSort);                                                              #PERCENT sort mutations according to radial SVM score (small to large)
my (@cnvUnsort, @cnvSort);                                                                          #PERCENT sort mutations according to cnv score (small to large)

my $givenPercentagePositionSVM;                                                                     #PERCENT length of all mutations with highest radial SVM score in each gene
my $givenPercentageRadialSVM;                                                                       #PERCENT radial SVM score at given percentile
my $givenPercentagePositionCNV;                                                                     #PERCENT length of all mutations with highest cnv score in each gene
my $givenPercentageCNV;                                                                             #PERCENT radial SVM score at given percentile



######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

&processArguments;
&formatConvert;
&loadDatabase;

######################################################################################################################################
#################################### seperately processing SNPs and CNVs #############################################################
######################################################################################################################################


open(OUT, "$outputFile") or die "iCAGES: cannot open input file $outputFile\n";

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



#################

    &processAnnovar;
    &processDrug;
    &getRelatedGenes;
    &generateOutput;

######################################################################################################################################
######################################################################################################################################
############################################################# subroutines ############################################################
######################################################################################################################################
######################################################################################################################################




######################################################################################################################################
###################################################### process arguments #############################################################
######################################################################################################################################



sub processArguments {
    GetOptions( 'help|h' => \$help,
    'percent|p=f' => \$percent,
    'manual|man|m' => \$manual,
    'downloadDB|d' => \$downloadDB,
    'id|i' => \$id
    )or pod2usage ();
    
    if($downloadDB){
        
        ################# commnad locations ######################
        $perlCommand = "$0";
        $perlCommand =~ /(.*)icages\.pl/;
        $icagesLocation = $1;
        $DBLocation = $icagesLocation. "db/";
        
        ######################## databases #######################
        $huiDB = $DBLocation . "phenolyzer.score";
        $cgcDB = $DBLocation . "cgc.gene";
        $keggDB = $DBLocation . "kegg.gene";
        $icagesDB = $DBLocation . "hg19_iCAGES.txt";
        $icagesIndex = $DBLocation . "hg19_iCAGES.idx";
        $funseqDB = $DBLocation . "hg19_funseq2.txt";
        $funseqIndex = $DBLocation . "hg19_funseq2.idx";
        $refGeneDB = $DBLocation . "hg19_refGene.txt";
        $refGeneIndex = $DBLocation . "hg19_refGeneMrna.fa";
        $cnvDB = $DBLocation . "hg19_cnv.txt";
        $suppressorDB = $DBLocation . "suppressor.gene";
        $oncogeneDB = $DBLocation . "oncogene.gene";
        $zscoreDB = $DBLocation . "drug.score";
        $biosystemDB = $DBLocation . "biosystem";
        
        unless(-e "$cgcDB"){
            !system("wget -O $cgcDB http://icages.usc.edu/download/icages/db/cgc.gene") or die "ERROR: cannot download database file $cgcDB\n";
        }
        
        unless(-e "$keggDB"){
            !system("wget -O $keggDB http://icages.usc.edu/download/icages/db/kegg.gene") or die "ERROR: cannot download database file $keggDB\n";
        }
        
        unless(-e "$huiDB"){
            !system("wget -O $huiDB http://icages.usc.edu/download/icages/db/phenolyzer.score") or die "ERROR: cannot download database file $huiDB\n";
        }
        
        unless(-e "$icagesDB"){
            !system("wget -O $icagesDB http://icages.usc.edu/download/icages/db/hg19_iCAGES.txt") or die "ERROR: cannot download database file $icagesDB\n";
        }
        
        unless(-e "$icagesIndex"){
            !system("wget -O $icagesIndex http://icages.usc.edu/download/icages/db/hg19_iCAGES.idx") or die "ERROR: cannot download database file $icagesIndex\n";
        }
        
        unless(-e "$refGeneDB"){
            !system("wget -O $refGeneDB http://icages.usc.edu/download/icages/db/hg19_refGene.txt") or die "ERROR: cannot download database file $refGeneDB\n";
        }
        
        unless(-e "$refGeneIndex"){
            !system("wget -O $refGeneIndex http://icages.usc.edu/download/icages/db/hg19_refGeneMrna.fa") or die "ERROR: cannot download database file $refGeneIndex\n";
        }
        
    }else{
        
        
        
        ######################## arguments ########################
        $percent=80 unless $percent;                                                                    #default value of percent
        $help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
        $manual and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
        @ARGV == 4 or @ARGV == 1 or pod2usage ();
        
        ################### configure file ########################
        
        $inputLocation = $ARGV[0];
        if(defined $ARGV[1]){
            $outputLocation = $ARGV[1];
            $outputLocation = $1 if $ARGV[1] =~ /(.*)\/$/;
        }else{
            $outputLocation = $1 if $ARGV[0] =~ /(.*)\/$/;
        }
        if(defined $ARGV[2]){
            $tempLocation = $ARGV[2];
            $tempLocation = $1 if $ARGV[2] =~ /(.*)\/$/;    
        }else{
            $tempLocation = $1 if $ARGV[0] =~ /(.*)\/$/;
        }
        if(defined $ARGV[3]){
            $logLocation = $ARGV[3];
                $logLocation = $1 if $ARGV[3] =~ /(.*)\/$/;    
        }else{
            $logLocation = $1 if $ARGV[0] =~ /(.*)\/$/;
        }
        $inputLocation = $1 if $ARGV[0] =~ /(.*)\/$/;
                
        ################# commnad locations ######################
        $perlCommand = "$0";
        $perlCommand =~ /(.*)icages\.pl/;
        $icagesLocation = $1;
        $DBLocation = $icagesLocation. "db/";
        
        #################### icages files #########################
        if(defined $id){
            $inputFile = $inputLocation . "/input-". $id . ".txt";
            $txt = $outputLocation . "/result-". $id . ".txt";
            $log = $logLocation . "/log-" . $id . ".log";
            $drug = $outputLocation . "/drug-". $id . ".txt";
            $suppressordrugs = $tempLocation . "/suppressordrugs-" . $id . ".drug";
            $oncogenedrugs = $tempLocation . "/oncogenedrugs-" . $id . ".drug";
            $otherdrugs = $tempLocation . "/otherdrugs-" . $id . ".drug";
            $relateddrugs = $tempLocation . "/relateddrugs-" . $id . ".drug";
            
        }else{
            $inputFile = $ARGV[0];
            $txt = $ARGV[0] . ".result.txt";
            $log = $ARGV[0] . ".log";
            $drug = $ARGV[0] . ".drug.txt";
            $suppressordrugs = $ARGV[0] . ".suppressordrugs.drug";
            $oncogenedrugs = $ARGV[0] . ".oncogenedrugs.drug";
            $otherdrugs = $ARGV[0] . ".otherdrugs.drug";
            $relateddrugs = $ARGV[0] . ".relateddrugs.drugs";
        }
        @drugs = ($suppressordrugs, $oncogenedrugs, $otherdrugs, $relateddrugs);
        
        ####################### ANNOVAR ##########################
        $annovarCommand = $icagesLocation . "bin/annovar";
        $dgidbCommand = $icagesLocation . "bin/DGIdb";
        $callDgidb = $dgidbCommand . "/getDrugList.pl";
        $callConvertToAnnovar = $annovarCommand . "/convert2annovar.pl";
        $callAnnovar = $annovarCommand . "/annotate_variation.pl";
        if(defined $id){
            $outputFile = $tempLocation ."/input-" .$id . ".annovar";                                   #define temp files for running ANNOVAR
        }else{
            $outputFile = $ARGV[0] . ".annovar";
        }
        $gene = $outputFile . ".variant_function";                                               #genes annotation from ANNOVAR
        $snpFile = $outputFile . ".snp";
        $cnvFile = $outputFile . ".cnv";
        $icages = $snpFile . ".hg19_iCAGES_dropped";                                                    #radial SVM score for each mutation
        #  $funseq = $snpFile . ".hg19_funseq2_dropped";
        $cnv = $cnvFile . ".hg19_cnv";
        $annovarLog = $outputFile .".log";
        
        
        ######################## databases #######################
        $huiDB = $DBLocation . "phenolyzer.score";
        $cgcDB = $DBLocation . "cgc.gene";
        $keggDB = $DBLocation . "kegg.gene";
        $icagesDB = $DBLocation . "hg19_iCAGES.txt";
        $icagesIndex = $DBLocation . "hg19_iCAGES.idx";
        # $funseqDB = $DBLocation . "hg19_funseq2.txt";
        # $funseqIndex = $DBLocation . "hg19_funseq2.idx";
        $refGeneDB = $DBLocation . "hg19_refGene.txt";
        $refGeneIndex = $DBLocation . "hg19_refGeneMrna.fa";
        $cnvDB = $DBLocation . "hg19_cnv.txt";
        $suppressorDB = $DBLocation . "suppressor.gene";
        $oncogeneDB = $DBLocation . "oncogene.gene";
        $zscoreDB = $DBLocation . "drug.score";
        $biosystemDB = $DBLocation . "biosystem";
        
        
        ######################## checking files #######################
        open(IN, "$inputFile") or die "iCAGES: cannot open input file $inputFile\n";
        open(SNP, ">$snpFile") or die "iCAGES: cannot open input file $snpFile\n";
        open(CNV, ">$cnvFile") or die "iCAGES: cannot open input file $cnvFile\n";
        open(CGC, "$cgcDB") or die "iCAGES: cannot open database file $cgcDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n";
        open(KEGG, "$keggDB") or die "iCAGES: cannot open database file $keggDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n";
        open(HUI, "$huiDB") or die "iCAGES: cannot open database file $huiDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n";
        open(SUP, "$suppressorDB") or die "iCAGES: cannot open database file $suppressorDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n";
        open(ONC, "$oncogeneDB") or die "iCAGES: cannot open database file $oncogeneDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n";
        open(ZSCORE, "$zscoreDB") or die "iCAGES: cannot open database file $zscoreDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n";
        open(BIO, "$biosystemDB") or die "iCAGES: cannot open database file $biosystemDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n";
        die "ERROR: cannot open database file $cnvDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n" unless -e $cnvDB ;
        die "ERROR: cannot open database file $icagesDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n" unless -e $icagesDB ;
        die "ERROR: cannot open index file $icagesIndex. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n" unless -e $icagesIndex ;
        die "ERROR: cannot open database file $refGeneDB. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n" unless -e $refGeneDB ;
        die "ERROR: cannot open database file $refGeneIndex. if you do not have this file in your .db/ folder, please download them using: icages.pl -downloadDB\n" unless -e $refGeneIndex ;
        
    }

}




######################################################################################################################################
######################################################## format convert ##############################################################
######################################################################################################################################


sub formatConvert{
    
    $formatCheckFirstLine = <IN>;
    
    print "NOTICE: start runing iCAGES packge at $nowString\n";
    print "NOTICE: start input file format checking and converting format if needed\n";
    push (@log, "iCAGES: start runing iCAGES packge at $nowString");
    push (@log, "iCAGES: start input file format checking and converting format if needed");
    
    chomp $formatCheckFirstLine;
    close IN;
    
    if($formatCheckFirstLine =~ /^##fileformat=VCF/){                                                               #VCF
        !system("$callConvertToAnnovar -format vcf4 $inputFile > $outputFile") or die "ERROR: cannot execute convert2annovar.pl for converting VCF file\n";
    }else{                                                                                                          #ANNOVAR
        !system("cp $inputFile $outputFile") or die "ERROR: cannot use input file $inputFile\n";
    }
}



######################################################################################################################################
######################################################## loading database ############################################################
######################################################################################################################################

sub loadDatabase{
    ###################### cgc ###########################
    push (@log, "iCAGES: start loading Cancer Gene Census data");
    print "NOTICE: start loading Cancer Gene Census data\n";
    
    while(<CGC>){
        chomp;
        my @line;
        @line = split;
        $cgc{$line[0]} = 1;
    }
    
    
    ###################### kegg ##########################
    push (@log, "iCAGES: start loading KEGG Cancer Pathway gene data");
    print "NOTICE: start loading KEGG Cancer Pathway gene data\n";
    
    while(<KEGG>){
        chomp;
        my @line;
        @line = split;
        $kegg{$line[0]} = $line[1];
    }
    
    
    ################### Phenolyzer #######################
    push (@log, "iCAGES: start extracting Phenolyzer score for each gene");
    print "NOTICE: start extracting Phenolyzer score for each gene\n";
    
    
    $minHui = 1;
    $maxHui = 0;
    
    while(<HUI>){
        chomp;
        my @line;
        @line = split(/\t/, $_);
        $hui{$line[0]} = $line[1];
        $maxHui = max($line[1], $maxHui);                                   #maxHui means the largest score in hui
        $minHui = min($line[1], $minHui);
    }
    
    
    
    ################### suppressor #######################
    push (@log, "iCAGES: start extracting suppressor genes");
    print "NOTICE: start extracting suppressor genes\n";
    
    
    while(<SUP>){
        chomp;
        $sup{$_} = 1;
    }
    
    close SUP;
    
    ################### oncogene #######################
    push (@log, "iCAGES: start extracting oncogenes");
    print "NOTICE: start extracting oncogenes\n";
    
    while(<ONC>){
        chomp;
        $onc{$_} = 1;
    }
    
    close ONC;
    
    ################### drug #######################
    push (@log, "iCAGES: start extracting drug scores");
    print "NOTICE: start extracting drug scores\n";
    
    
    
    while(<ZSCORE>){
        chomp;
        my @line;
        @line = split(/\t/, $_);
        $zscore{$line[0]} = $line[1];
    }
    close ZSCORE;
    
}



######################################################################################################################################
################################################ processing ANNVOAR output ###########################################################
######################################################################################################################################

sub processAnnovar{
    
    open(TXT, ">$txt") or die "ERROR: cannot open file $txt\n";
    open(LOG, ">$log") or die "ERROR: cannot open file $log\n";
    open(DRUG, ">$drug") or die "ERROR: cannot open file $drug\n";
    open(GENE, "$gene") or die "ERROR: cannot open file $gene\n";
    open(ANNLOG, "$annovarLog") or die "ERROR: cannot open file $annovarLog\n";
    open(MYCNV, "$cnv") or die "ERROR: cannot open file $cnv\n";
    open(RADIAL, "$icages") or die "ERROR: cannot open file $icages\n";
    # open(FUNSEQ, "$funseq") or die "ERROR: cannot open file $funseq\n";
    
    
    ############# log and print to iCAGES ##############
    while(<ANNLOG>){
        chomp;
        push (@log, $_);
    }
    
    for(0..$#log){
        print LOG "$log[$_]\n";
    }
    
    
    ############### radial SVM score ###################
    
    print LOG "iCAGES: start extracting radial SVM score for each mutation from ANNOVAR OUTPUT\n";
    print "NOTICE: start extracting radial SVM score for each mutation from ANNOVAR OUTPUT\n";
    
    while(<RADIAL>){
        chomp;
        my @line;
        my $key;
        @line = split(/\t/, $_);
        $key = "$line[2]:$line[3]:$line[4]:$line[5]:$line[6]";
        $radialSVM{$key} = $line[1];
    }
    
    
    ############### funseq score ###################
    
    # print LOG "iCAGES: start extracting funseq score for each mutation from ANNOVAR OUTPUT\n";
    # print "NOTICE: start extracting funseq score for each mutation from ANNOVAR OUTPUT\n";
    
    # while(<FUNSEQ>){
    #    chomp;
    #    my @line;
    #    my $key;
    #    @line = split(/\t/, $_);
    #    $key = "$line[2]:$line[3]:$line[4]:$line[5]:$line[6]";
    #    $funseq{$key} = $line[1];
    # }
    
    
    
    
    ############# structural variation ################
    
    print LOG "iCAGES: start processing structural variations from ANNOVAR OUTPUT\n";
    print "NOTICE: start processing structural variations from ANNOVAR OUTPUT\n";
    
    while(<MYCNV>){
        chomp;
        my @line;
        my $key;
        my $score;
        @line = split(/\t/, $_);
        $score = $line[1];
        $score =~ /Score=(.*);/;
        $score = $1;
        $key = "$line[2]:$line[3]:$line[4]:$line[5]:$line[6]";
        $cnv{$key} = $score;
    }
    
    
    
    
    ############### gene annotation ###################
    
    print LOG "iCAGES: start parsing gene annotation files generated from ANNOVAR output\n";
    print "NOTICE: start parsing gene annotation files generated from ANNOVAR output\n";
    
    
    
    while(<GENE>){
        chomp;
        my @line;
        my $key;                                                            #hash key used for fetch radial SVM score from %radialSVM: mutation->radialSVM
        my $gene;
        @line = split(/\t/, $_);
        $key = "$line[2]:$line[3]:$line[4]:$line[5]:$line[6]";
        $gene = $line[1];
        next unless defined $gene;
        
        if(exists $radialSVM{$key}){
            if( exists $icages{$gene}{"radialSVM"}){
                $icages{$gene}{"radialSVM"} = max($icages{$gene}{"radialSVM"}, $radialSVM{$key});
            }else{
                $icages{$gene}{"radialSVM"} = $radialSVM{$key};
            }
        }else{
            $icages{$gene}{"radialSVM"} = 0;
        }
        
        if(exists $cnv{$key} and exists $sup{$gene}){
            if( exists $icages{$gene}{"cnv"}){
                $icages{$gene}{"cnv"} = max($icages{$gene}{"cnv"}, $cnv{$key});
            }else{
                $icages{$gene}{"cnv"} = $cnv{$key};
            }
        }else{
            $icages{$gene}{"cnv"} = 0;
        }
        
        # if(exists $funseq{$key}){
        #    if( exists $icages{$gene}{"funseq"}){
        #        $icages{$gene}{"funseq"} = max($icages{$gene}{"funseq"}, $radialSVM{$key});
        #    }else{
        #        $icages{$gene}{"funseq"} = $funseq{$key};
        #    }
        #}
        print "$gene:$hui{$gene}:\n";
        if(exists $hui{$gene}){
            if( exists $icages{$gene}{"phenolyzer"}){
                $icages{$gene}{"phenolyzer"} =  max( $icages{$gene}{"phenolyzer"}, $hui{$gene});
            }else{
                $icages{$gene}{"phenolyzer"} =  $hui{$gene};
            }
        }else{
            $icages{$gene}{"phenolyzer"} = 0;
        }
    }
    
    foreach my $gene (sort keys %icages){
        $driver{$gene} = 1;
    }
    close GENE;
}

sub getRelatedGenes{
    print LOG "iCAGES: start getting list of genes related to driver genes\n";
    print "NOTICE: start getting list of genes related to driver genes\n";
    
    while(<BIO>){
        chomp;
        my @line = split(/\t/, $_);
        next unless exists $line[0];
        next unless exists $line[1];
        next unless exists $line[2];
        if(exists $driver{$line[0]}){
            $related{$line[1]}{$line[0]} = $line[2];
            push @relatedgenes, $line[1];
        }
    }
    close BIO;
}



######################################################################################################################################
################################################# processing drugs output ############################################################
######################################################################################################################################

sub processDrug{
    
    ############### get drugs according to the list of cancer driver genes ###################
    
    print LOG "iCAGES: start getting list of drugs for predicted cancer driver genes\n";
    print "NOTICE: start getting list of drugs for predicted cancer driver genes\n";
    
    my ($oncogenes, $suppressors, $othergenes, $relatedgenes);


    $suppressors = join(",", @suppressors);
    $oncogenes = join(",", @oncogenes);
    $othergenes = join(",", @othergenes);
    $relatedgenes = join(",", @relatedgenes);
    
    if($suppressors ne ""){
        !system("$callDgidb --genes='$suppressors' --interaction_type='activator,other/unknown,n/a,inducer,stimulator' --source_trust_levels='Expert curated' --output='$suppressordrugs'") or die "ERROR: cannot get drugs\n";
    }
    if($oncogenes ne ""){
        !system("$callDgidb --genes='$oncogenes' --interaction_type='inhibitor,suppressor,antibody,antagonist,blocker,other/unknown,n/a' --source_trust_levels='Expert curated' --output='$oncogenedrugs'") or die "ERROR: cannot get drugs\n";
    }
    if($othergenes ne ""){
        !system("$callDgidb --genes='$othergenes' --interaction_type='inhibitor,suppressor,antibody,antagonist,blocker,activator,other/unknown,n/a,inducer,stimulator' --source_trust_levels='Expert curated' --output='$otherdrugs'") or die "ERROR: cannot get drugs\n";
    }
    if($relatedgenes ne ""){
        !system("$callDgidb --genes='$relatedgenes' --interaction_type='inhibitor,suppressor,antibody,antagonist,blocker,activator,other/unknown,n/a,inducer,stimulator' --source_trust_levels='Expert curated' --output='$relateddrugs'") or die "ERROR: cannot get drugs\n";
    }
    
    ######################### get drugs zscores #######################

    
    for(0..$#drugs-1){
        open(DRUGOUT, "$drugs[$_]") or next;
        while(<DRUGOUT>){
            chomp;
            my @line = split(/\t/, $_);
            my $drugscore ;
            if(exists $zscore{$line[1]}){
                $drugscore = $zscore{$line[1]};
                $drug{$line[0]}{$line[1]} = $drugscore;
            }else{
                $drug{$line[0]}{$line[1]} = 0;
            }
        }
    }
   
    foreach my $key (sort keys %drug){
        foreach my $drugkey (sort keys %{$drug{$key}}){
            print DRUG "$key\t$drugkey\t$drug{$key}{$drugkey}\n";
            
        }
    }
    
    close DRUGOUT;
    close DRUG;
}



######################################################################################################################################
################################################# generating final output ############################################################
######################################################################################################################################

sub generateOutput{
    
    print LOG "iCAGES: start generating output\n";
    print "NOTICE: start generating output\n";
    
    $iGene = 0;
    
    
    # print TXT "GeneName\tRadialSVMScore\tCNV\tFunseq\tPhenolyzerSocre\n";
     print TXT "GeneName\tRadialSVMScore\tCNV\tPhenolyzerSocre\n";
    
    foreach my $gene (sort keys %icages){              #loop to get information for each gene
        #   foreach my $key ("radialSVM", "cnv", "funseq", "phenolyzer"){
            $icages{$gene}{"radialSVM"} = "\"NA\"" unless exists $icages{$gene}{"radialSVM"};
            $icages{$gene}{"cnv"} = "\"NA\"" unless exists $icages{$gene}{"cnv"};
            $icages{$gene}{"phenolyzer"} = "\"NA\"" unless exists $icages{$gene}{"phenolyzer"};
            # print TXT "$gene\t$icages{$gene}{\"radialSVM\"}\t$icages{$gene}{\"cnv\"}\t$icages{$gene}{\"funseq\"}\t$icages{$gene}{\"phenolyzer\"}\n";
            print TXT "$gene\t$icages{$gene}{\"radialSVM\"}\t$icages{$gene}{\"cnv\"}\t$icages{$gene}{\"phenolyzer\"}\n";
    };
    
    $nowString = localtime;

    print LOG "iCAGES: finished at $nowString\n";
    print "NOTICE: finished at $nowString\n";
    
}



######################################################################################################################################
############################################################ manual page #############################################################
######################################################################################################################################

=head1 NAME                                                                                                                                                                                                                                                                 
           
 iCAGES (integrated CAncer GEnome Score) command line package for web interface.
                                                              
=head1 SYNOPSIS
                                                                                                                                                                                                      
 icages.pl [options] <input>                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                               
 Options:                                                                                                                                                                                                                                                                     
        -h, --help                      print help message   
        -m, --manual                    print manual message
        -p, --percent <float>           percentile above which deleterious mutations are selected to determine the candidate cancer driver genes for each patient (default: 80)          
        

 Function: iCAGES predicts cancer driver genes given protein altering somatic mutations from a patient.                                                                                                                                                                                                                                                                                                                                                                                
 Example: icages.pl /path/to/input.txt 
                                                                                                                                                                                                                                                                               
 Version: 1.0

 Last update: Wed Aug 20 19:47:38 PDT 2014
 
=head1 OPTIONS

=over 8

=item B<--help>

 print a brief usage message and detailed explanation of options.

=item B<--manual>

 print the manual page and exit.

=back

=cut

