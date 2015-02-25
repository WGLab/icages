#!/usr/bin/perl
use strict;
use warnings;

######################################################################################################################################
######################################################## variable declaration ########################################################
######################################################################################################################################

my $icagesLocation = "$0";
$icagesLocation =~ /(.*)icages\.pl/;
$icagesLocation = $1;
$icagesLocation = "./" if $icagesLocation eq "";
$icagesLocation =~ /(.*)\/$/;
$icagesLocation = $1;


######################################################################################################################################
########################################################### main  ####################################################################
######################################################################################################################################

print "NOTICE: start initiating iCAGES package. This process may take a long time and requires more than 200G of storage.\n";

############################ download database file ############################
!system("wget -O $icagesLocation/db.tar.gz http://icages.usc.edu/download/icages/db.tar.gz") or die "NOTICE: cannot download database file\n";

########################### uncompress database file ###########################
!system("tar -zxvf $icagesLocation/db.tar.gz") or die "NOTICE: cannot uncompress database file\n";

############################# download script file #############################
!system("mkdir $icagesLocation/bin/DGIdb");
!system("mkdir $icagesLocation/bin/annovar");
!system("wget -O $icagesLocation/bin/annovar/convert2annovar.pl http://icages.usc.edu/download/icages/bin/annovar/convert2annovar.pl") or die "NOTICE: cannot download script file\n";
!system("wget -O $icagesLocation/bin/annovar/annotate_variation.pl http://icages.usc.edu/download/icages/bin/annovar/annotate_variation.pl") or die "NOTICE: cannot download script file\n";
!system("wget -O $icagesLocation/bin/DGIdb/getDrugList.pl http://icages.usc.edu/download/icages/bin/DGIdb/getDrugList.pl") or die "NOTICE: cannot download script file\n";
!system("wget -O $icagesLocation/bin/DGIdb/JSON.pm http://icages.usc.edu/download/icages/bin/DGIdb/JSON.pm") or die "NOTICE: cannot download script file\n";


############################# install package file #############################
!system("mkdir $icagesLocation/bin/myperllib");

# cpanm
if(!testFunction($icagesLocation, "cpanm") and !testFunction($icagesLocation, "~/perl5/bin/cpanm")){
    if(-e "~/perl5"){
        !system("rm -fr ~/perl5");
    };
    !system("curl -L http://cpanmin.us | perl - App::cpanminus");
    writeToBash("export PATH=$PATH:.:$HOME/perl5/bin");
}

# local::lib
!system("cpanm --local-lib=$icagesLocation/bin/myperllib local::lib");
writeToBash("eval $(perl -I \"$icagesLocation/bin/myperllib/lib/perl5\" -Mlocal::lib)");

# JSON, HTTP::Request, LWP
!system("cpanm JSON");
!system("cpanm HTTP::Request");
!system("cpanm LWP");


######################################################################################################################################
######################################################## subroutines  ################################################################
######################################################################################################################################


sub testFunction(){
    my $icagesLocation = shift;
    my $commnd = shift;
    my $testingFile = "$icagesLocation/comman_testing.txt";
    my $testFileContent;
    !system("$commnd && echo OK || echo Failed > $testingFile")
    open(TEST, $testingFile) or die;
    $testFileContent = <TEST>;
    chomp $testingFileContent;
    if($testingFileContent eq "OK"){
        return 1;
    }else{
        return 0;
    };
}



sub writeToBash(){
    my $content = shift;
    my $bashfile;
    if(-e "~/.bash_profile"){
        $bashfile = "~/.bash_profile";
    }else(-e "~/.bashrc"){
        $bashfile = "~/.bashrc";
    }
    
    if(-e $bashfile){
        !system("echo \"$content\" >> $bashfile");
    }else{
        !system("touch $bashfile");
        !system("echo '# .bashrc' >> $bashfile");
        !system("echo \"$content\" >> $bashfile");
    }
    !system("source $bashfile");
}



