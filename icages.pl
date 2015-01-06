#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################
my ($inputLocation, $icagesLocation);
my ($icagesMutation, $icagesGene, $icagesDrug);

######################################################################################################################################
############################################################# main  ##################################################################
######################################################################################################################################
($inputLocation, $icagesLocation) =  &processArguments();
$icagesMutation = "$inputLocation/bin/icagesMutation.pl";
$icagesGene = "$inputLocation/bin/icagesGene.pl";
$icagesDrug = "$inputLocation/bin/icagesDrug.pl";

!system("perl $icagesMutation $inputLocation $icagesLocation") or die "ERROR: cannot call icagesMutation module\n";
!system("perl $icagesGene $inputLocation $icagesLocation") or die "ERROR: cannot call icagesGene module\n";
!system("perl $icagesDrug $inputLocation $icagesLocation") or die "ERROR: cannot call icagesDrug module\n";

######################################################################################################################################
########################################################## subroutines ###############################################################
######################################################################################################################################

sub processArguments {
    my ($help, $manual, $inputLocation, $icagesLocation);
    GetOptions( 'help|h' => \$help,
    'manual|man|m' => \$manual
    )or pod2usage ();
    ######################## arguments ########################
    $help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
    $manual and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
    @ARGV == 1 or pod2usage ();
    ################### locations ########################
    $inputLocation = $1 if $ARGV[0] =~ /(.*)\/$/;
    $icagesLocation = "$0";
    $icagesLocation =~ /(.*)icages\.pl/;
    $icagesLocation = $1;
    return ($inputLocation, $icagesLocation);
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
 
 Function: iCAGES predicts cancer driver genes given somatic mutations (in ANNOVAR/VCF format) from a patient.
 Example: icages.pl /path/to/input.vcf
                                                                                                                                                                                                                                                                               
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

