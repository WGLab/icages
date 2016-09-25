#!/usr/bin/perl
use strict;
use warnings;
use Fcntl qw/SEEK_END SEEK_SET/;
use Data::Dumper;
use Time::HiRes qw/gettimeofday tv_interval/;

my $BINSIZE = 10; #bin size for indexing
my $ENDOFFILE = "END_OF_GENOMELOCUSFINDER_INDEX";
my $DEBUG = 0;
my $suffix = "binsize${BINSIZE}gidx";
my $usage = "Usage: $0 
index	index a sorted tab-delmited text file
find 	retrieve a record for some loci\n";
die $usage unless @ARGV >= 1;
my $subprogram = shift @ARGV;
if ($subprogram eq 'index') {
    &index(@ARGV);
} elsif ($subprogram eq 'find') {
    &find(@ARGV);
} else {
    die $usage;
}

##########WARNING##################################################################
#typical use cases for this design
#small vocabulary in 1st column
#large database file
#large query file
#only a few annotations for each genomic locus
#not many overlappic loci (e.g. no chr1:1-10, chr1:1-100, chr1:1-1000 etc)
####MIGHT NOT WORK WELL UNDER OTHER CIRCUMSTATNCES################################
##################SUBROUTINES#######################
sub index {
    #generate hierarchical index for a sorted tab-delimited txt file
    my $usage = "$0 index <fasta index file> <sorted BED file>\n";
    die $usage unless @_ == 2;
    warn "Start indexing...\n";
    my $t0 = [gettimeofday];
    my $fai = shift;
    my $in = shift;
    my $out = "$in.$suffix";
    #input format
    #1	10001	10001	T	A	0.18521432
    #types of binary representation: Q: unsigned 64-bit integer, 8 bytes
    #index layout
    #############
    #bin index: 8-byte unsigned integer marking start of each $BINSIZE bp bin in the file being indexed
    #key index: (number of keys)*8 bytes of unsigned integers marking start position of each bin index 
    #		section and length associated with each key (for chromosomes, this is equivalent to chr
    #		length), the order of key indeces is same as the order of vocabulary
    #vocabulary: a list of keys separated by tabs 
    #offset: 8 bytes unsiged integer, mark total bytes for vocabulary
    #######################An example index file################################
    #-----------------------------------------------------------------
    #(bin index section)
    #(chr1 bins)
    #0x0000 0000 0000 0001
    #0x0000 0000 0000 00FA
    #...
    #(chr2 bins)
    #...
    #-----------------------------------------------------------------
    #(key index section)
    #(chr1 bin start in this index)0x00000001 ||(chr1 length)0xE9A6740
    #(chr2 bin start in this index)0x000000EB ||(chr2 length)0xE7EED8D
    #...
    #-----------------------------------------------------------------
    #(vocabulary section)
    #chr1\tchr2\t...
    #-----------------------------------------------------------------
    #(the offset)0x000010AE
    #-----------------------------------------------------------------
    my %key_length = &readFastaIndex($fai);
    my $vocabulary_string;
    my %key_start;
    my $offset;

    open IN,'<',$in or die "$in: $!";
    open OUT,'>',$out or die "$out: $!";
    binmode OUT; #index is in binary format for space efficiency

    my $current_db_position = 0; #current position in database file
    my $count_total = 0;
    my $previous_key;
    my ($previous_start, $previous_end);
    my $walker = 0; #measures how far we have walked on each chr, in unit of $BINSIZE
    #every time we pass $BINSIZE, record position
    #reset for each new chr
    while (<IN>) {
	$count_total++;
	warn "NOTICE: $count_total records processed in ".tv_interval($t0)." seconds\n" if $count_total % 1_000_000 == 0;
	#example database record
	#0	1	2	3	4	5
	#1	10001	10001	T	A	0.18521432
	#assume database is sorted
	chomp;
	my @f = split /\t/;
	die "ERROR: expect at least 5 fields at line $. of $in\n" unless @f >= 5;
	die "ERROR: Database $in not sorted in ascending order at line $. of $in.\n" if defined $previous_start and defined $previous_key and
	$previous_key eq $f[0] and (($f[1] < $previous_start) or ($previous_start == $f[1] and $f[2] < $previous_end));
	die "ERROR: start larger than end at line $. of $in\n" unless $f[2] >= $f[1];
	die "ERROR: range out of predefined bounds of $f[0] at line $. of $in.\n"
	if $key_length{$f[0]} < $f[1] or $key_length{$f[0]} < $f[2];

	if (!defined $previous_key or $previous_key ne $f[0]) {
	    #we got a new key
	    $key_start{$f[0]} = tell OUT;
	    if (defined $previous_key) {
		#modify length according to how far we have walked
		$key_length{$previous_key} = $walker;
	    }
	    $walker = 0;
	}

	while ($f[1] >= $walker) {
	    print $current_db_position,"\n" if $DEBUG;
	    print OUT pack("Q",$current_db_position);
	    $walker += $BINSIZE;
	}

	$current_db_position = tell IN;
	$previous_key = $f[0];
	($previous_start,$previous_end) = @f[1,2];
    }
    close IN;
    #record last key's end position
    $key_length{$previous_key} = $walker;

    #print out key starts and lengths
    for my $key(sort keys %key_start) {
	print OUT pack("Q",$key_start{$key});
	print OUT pack("Q",$key_length{$key});
	print $key_start{$key},"\n" if $DEBUG;
	print $key_length{$key},"\n" if $DEBUG;
    }
    #print vocabulary
    $vocabulary_string = join("\t",sort keys %key_start);
    print OUT pack("A".length($vocabulary_string), $vocabulary_string);
    print $vocabulary_string,"\n" if $DEBUG;

    #print offset
    $offset = length($vocabulary_string);
    print OUT pack("Q", $offset);
    print $offset,"\n" if $DEBUG;

    #print EOF mark
    print OUT pack("A".length($ENDOFFILE), $ENDOFFILE);
    close OUT;

    warn "Indexing done. Bin size: $BINSIZE. Assume ASCII Encoding for characters.\n";
    warn "Processed $count_total records in ".tv_interval($t0)." seconds\n";
}
sub find {
    #retrieve recording overlapping with queries
    my $usage = "$0 find <database> <query file> <database name> <output file>\n";
    die $usage unless @_ == 4;
    my $t0 = [gettimeofday];
    my $db = shift;
    my $index_file = "$db.$suffix";
    my $query = shift;
    my $databaseName = shift;
    my $outputFile = shift;
    warn "Start querying $db...\n";
    #how to find the record for a specific locus?
    #first read the offset at the last 4 bytes
    #then load the vocabulary into hash in memory (size is determined by offset)
    #then locate key index (each key has one start 4 bytes, one length 4 bytes, total 8 bytes).
    #store the key indeces in hash
    #locate bin index by key index and genomic coordinate
    #go to the file, read from bin start until no more records
    my %key_index = &loadKeyIndex($index_file);
    #example query
    #1	114438528	114438528	G	A	TCGA-06-0155-01B-01D-1492-08
    open DB,'<',$db or die "$db: $!";
    open INDEX,'<',$index_file or die "$index_file: $!";
    open QUERY,'<',$query or die "$query: $!";
    open OUT, ">", $outputFile or die "$outputFile: $!";
    
    my $count_match = 0;
    my $count_total = 0;
    my $count_total_db_lookup = 1;

    while (<QUERY>) {
	$count_total++;
	warn "NOTICE: $count_total records processed in ".tv_interval($t0)." seconds\n" if $count_total % 1_000_000 == 0;
#	print "query: $_" ;
	chomp;
	my @f = split /\t/;
	die "ERROR: expect at least 5 fields at line $. of $query\n" unless @f >= 5;
	warn "QUERY: @f\n" if $DEBUG;
	my ($chr, $start, $end, $ref, $alt) = @f[0..4];
	# chr start
	if ($chr =~ /chr/) {
	    $chr =~ /chr(.*)/;
	    $chr = $1;
	}
	if ((not exists $key_index{$chr}) or $start > $key_index{$chr}->[1]) {
	    warn "NOTICE: @f no match found\n" if $DEBUG;
	} else {
	    seek INDEX, $key_index{$chr}->[0] + 8*(int $start/$BINSIZE), SEEK_SET;
	    my $position_in_db;
	    read INDEX, $position_in_db, 8;
	    $position_in_db = unpack("Q", $position_in_db);
	    seek DB, $position_in_db, SEEK_SET;
	    while (<DB>) {
		#search within database
		#until no more annotations can be possibly found
		$count_total_db_lookup++;
#		print;
		chomp;
		my @db_fields = split /\t/;
		warn "DB:@db_fields\n" if $DEBUG;


		if ($db_fields[0] eq $chr and $db_fields[1] == $start and $db_fields[2] == $end
			and $db_fields[3] eq $ref and $db_fields[4] eq $alt) {
		    print OUT join("\t", $databaseName, $db_fields[5], @f[0..4]), "\n";
		    $count_match++;
		}
		if ($db_fields[0] ne $chr or $db_fields[1] > $start) {
		    last; #if chromosomes are not equal, or the start coordinate is larger than query start
		    #there is no need to continue searching
		}
	    }
	}
    }
    close QUERY;
    close DB;
    close INDEX;
    warn "Querying done. Found $count_match matches in $count_total queries.\n";
    warn "Processed $count_total records in ".tv_interval($t0)." seconds\n";
    warn "Average database lookup: ".($count_total_db_lookup/$count_total)."\n";
}
sub loadKeyIndex {
    #load vocabulary and key indeces of the index file
    my $idx = shift;

    my $offset;
    my %key_index; #values are [start_position, length]
    my @vocabulary;
    my $buffer;
    open IN,'<',$idx or die "$!";
    binmode IN;

    #check end
    seek IN, -length($ENDOFFILE), SEEK_END or die "$!";
    read IN, $buffer, length($ENDOFFILE) or die "$!";
    $buffer = unpack("A".length($ENDOFFILE), $buffer);
    die "Incomplete or corrupted index file, please re-index the database.\n" unless $buffer eq $ENDOFFILE;

    #offset
    seek IN, -8-length($ENDOFFILE), SEEK_END or die "$!";
    read IN, $offset, 8 or die "$!";
    $offset = unpack("Q",$offset);
    print "offset: $offset\n" if $DEBUG;

    #vocabulary
    ##reopen the filehandle because last time we read to the end
    #open IN,'<',$idx or die "$!";
    #binmode IN;
    seek IN, -(8 + $offset + length($ENDOFFILE)), SEEK_END or die "$!";
    read IN, $buffer, $offset;
    @vocabulary = split /\t/,(unpack ("A".$offset,$buffer));
    print "@vocabulary" if $DEBUG;

    #key index
    seek IN, -(length($ENDOFFILE) + 8 + $offset + 16*@vocabulary), SEEK_END;
    for my $key(@vocabulary) {
	my ($start, $len);
	read IN, $start, 8;
	read IN, $len, 8;
	$start = unpack("Q", $start);
	$len = unpack("Q", $len);
	$key_index{$key} = [$start, $len];
    }

    print Dumper(%key_index) if $DEBUG;
    close IN;
    return %key_index;
}
sub readFastaIndex {
    my $in = shift;
    #chr1	243000000	80	120
    my %vocabulary;
    open IN,'<',$in or die "$in: $!";
    while (<IN>) {
	my @f = split /\t/;
	$vocabulary{$f[0]} = $f[1];
    }
    close IN;
    return %vocabulary;
}
