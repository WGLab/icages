#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################
my ($icagesMutation, $icagesGene, $icagesDrug, $icagesJson);
my ($inputFile, $inputDir, $icagesLocation, $tumor, $germline, $id, $subtype , $logDir, $outputDir, $tempDir, $prefix, $bed, $hg19);

######################################################################################################################################
############################################################# main  ##################################################################
######################################################################################################################################
($inputDir, $icagesLocation, $tumor, $germline, $id, $subtype , $logDir, $outputDir, $tempDir, $prefix, $bed, $hg19) =  &processArguments();
&checkReady($icagesLocation);
$icagesMutation = $icagesLocation. "bin/icagesMutation.pl";
$icagesGene = $icagesLocation . "bin/icagesGene.pl";
$icagesDrug = $icagesLocation . "bin/icagesDrug.pl";
$icagesJson = $icagesLocation . "bin/icagesJson.pl";
$inputFile = $ARGV[0];
!system("perl $icagesMutation $inputFile $inputDir $icagesLocation $tumor $germline $id $prefix $bed $hg19") or die "ERROR: cannot call icagesMutation module\n";
!system("perl $icagesGene $inputDir $icagesLocation $subtype $prefix ") or die "ERROR: cannot call icagesGene module\n";
!system("perl $icagesDrug $inputDir $icagesLocation $prefix") or die "ERROR: cannot call icagesDrug module\n";
!system("perl $icagesJson $inputDir $icagesLocation $prefix") or die "ERROR: cannot call icagesJson module\n";
&moveFiles($inputDir, $prefix, $logDir, $outputDir, $tempDir);

######################################################################################################################################
########################################################## subroutines ###############################################################
######################################################################################################################################


sub moveFiles{
    my ($inputDir, $prefix, $logDir, $outputDir, $tempDir);
    my ( $outputFile, $tempFile, $logFile, $jsonFile);
    $inputDir = shift;
    $prefix = shift;
    $logDir = shift;
    $outputDir = shift;
    $tempDir = shift;
    if($inputDir ne $outputDir){
	$outputFile = "$inputDir/$prefix*icages*.csv";
	$jsonFile = "$inputDir/$prefix*.json";
	!system("mv $outputFile $outputDir") or die "ERROR: cannot move iCAGES file\n";
	!system("mv $jsonFile $outputDir") or die "ERROR: cannot move iCAGES file\n";
    }
    if($inputDir ne $logDir){
	$logFile = "$inputDir/$prefix*.log";
	!system("mv $logFile $logDir") or die "ERROR: cannot move iCAGES file\n";	
    }
    if($inputDir ne $tempDir){
	$tempFile = "$inputDir/$prefix.*";
	!system("mv $tempFile $tempDir") or die "ERROR: cannot move iCAGES file\n";
    }
}

sub checkReady() {
    my $icagesLocation = shift;
    my $dbLocation = $icagesLocation . "db";
    if(-d $dbLocation){
        return 1;
    }else{
        die "ERROR: please download iCAGES database first https://github.com/WangGenomicsLab/icages \n";
    }
}

sub processArguments {
    my ($help, $manual, $tumor, $germline, $id, $subtype, $logDir, $outputDir, $tempDir, $prefix, $inputDir, $inputLocation, $icagesLocation, $bed, $hg);
    ################### initialize arguments ##################                                                                                                                                           
    GetOptions( 'help|h' => \$help,
    'manual|man|m' => \$manual,
    'tumor|t=s' => \$tumor,   # name for tumor in the vcf file                                                                                                                                             
    'germline|g=s' => \$germline ,   # name for germline in the vcf file                                                                                                                                  
    'id|i=s' => \$id,   # sample identifier for the person of interest for multiple sample vcf file                                                                                                        
    'subtype|s=s' => \$subtype, # cancer subtype                                                                                                                                        
    'logdir=s' => \$logDir, # log directory                                                                                                                                                                
    'outputdir=s' => \$outputDir,
    'tempdir=s' => \$tempDir,
    'prefix|p=s' => \$prefix,
    'bed|b=s' => \$bed, # bed file describing structural variations
    'buildver=s' => \$hg		
	)or pod2usage ();
    ################### locations ########################
    if($hg and $hg ne "hg19" and $hg ne "hg38" and $hg ne "hg18"){
	pod2usage ();
    }
    @ARGV == 1 or pod2usage (); # check only has one argument 
    $inputLocation = $ARGV[0];
    $inputDir = $inputLocation;
    if($inputDir =~ /\//){
	$inputDir =~ /(.*\/)(.*?)$/;
	$inputDir = $1;
    }else{
	$inputDir = "./" ;
    }
    $icagesLocation = "$0";
    $icagesLocation =~ /(.*)icages\.pl/;
    $icagesLocation = $1;
    $icagesLocation = "./" if $icagesLocation eq "";
    ###### all directories should end up with / ###
    if(!$prefix){
	$prefix = $ARGV[0];
	if($prefix =~ /\//){
	    $prefix =~ /(.*\/)(.*?)$/;   
	    $prefix = $2;
	}
    }
    if(!$tumor){
	$tumor = "NA";
    }
    if(!$germline){
	$germline = "NA";
    }
    if(!$id){
	$id = "NA";
    }
    if(!$subtype){
	$subtype = "NA";
    }
    if(!$bed){
	$bed = "NA";
    }
    if(!$hg){
	$hg = "hg19";
    }
    if(!$logDir){
	if(!$outputDir){
	    $logDir = $inputDir;
	}else{
	    $logDir = $outputDir;
	}
    }
    if(!$outputDir){
	$outputDir = $inputDir;
    }
    if(!$tempDir){
	if(!$outputDir){
	    $tempDir = $inputDir;
	}else{
	    $tempDir = $outputDir;
	}
    }
    if(-d $logDir){
        if(!($logDir =~ /\/$/)){
	    $logDir = $logDir . "/";
	}	
    }else{
	if(!($logDir =~ /\/$/)){
            $logDir = $logDir . "/";
        }
	mkdir($logDir) or die "ERROR: no such directory for log files\n";
    }
    if(-d $outputDir){
	if(!($outputDir =~ /\/$/)){
            $outputDir = $outputDir ."/";
        }
    }else{
	if(!($outputDir =~ /\/$/)){
            $outputDir = $outputDir ."/";
        }
	mkdir($outputDir) or die "ERROR: no such directory for output files\n";
    }
    if(-d $tempDir){
	if(!($tempDir =~ /\/$/)){
            $tempDir = $tempDir ."/";
        }
    }else{
	if(!($tempDir =~ /\/$/)){
            $tempDir = $tempDir ."/";
        }
	mkdir($tempDir) or die "ERROR: no such directory for temp files\n";
    }

    ######################## arguments ########################
    $help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
    $manual and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
    return ($inputDir, $icagesLocation, $tumor, $germline, $id, $subtype , $logDir, $outputDir, $tempDir, $prefix, $bed, $hg);
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
        -t, --tumor                     name of column that contains tumor mutations in your vcf file (if you have multiple samples with tumor mutations, please use this option to select tumor mutations that you want to analyze)
        -g, --germline                  name of column that contains germline mutations in your vcf file (if you have multiple samples with germline mutations, please use this option to select germline mutations that you want to compare your tumor mutations against to generate somatic mutation profiles for the sample you want to analyze)
        -i, --id                        name of column that contains somatic mutations in your multiple sample vcf file with only somatic mutations (if you have multiple samples with tumor and germline mutations, please use -g and -t options instead)
        -s, --subtype                   subtype of the cancer, valid options include "ACC", "BLCA", "BRCA", "CESC", "CHOL", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "STAD", "TGCT", "TGCA", "THYM", "UCEC", "UCS", "UVM"
        --logdir                        directory for log files generated by iCAGES
        --tempdir                       directory for temporary files generated by iCAGES
        --outputdir                     directory for output files generated by iCAGES
        -p, --prefix                    prefix of all files generated by iCAGES
        -b, --bed                       additional bed file specifying the location of structural variations in the sample
        --buildver                      reference genome version, valid options include "hg19" (default), "hg38" and "hg18"
 
 Function: iCAGES predicts cancer driver genes given somatic mutations (in ANNOVAR/VCF format) from a patient.
 
 Example: icages.pl /path/to/input.vcf
 
 Installation: before using iCAGES, please first install it by 'perl icagesInitiate.pl' command.
 
 Version: 1.0
 
 Last update: Wed Feb 25 12:51:17 PST 2015
 
=head1 OPTIONS

=over 8

=item B<--help>

 print a brief usage message and detailed explanation of options.

=item B<--manual>

 print the manual page and exit.

=back

=cut

