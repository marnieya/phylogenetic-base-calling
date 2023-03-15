#!/usr/bin/perl
use strict;
use warnings;
if ( scalar(@ARGV) < 6){
	print "./simulate_paired_end_reads.pl\n";
	print "ARGV[0]: Fasta to simulate reads from\n";
	print "ARGV[1]: Prefix\n";
	print "ARGV[2]: Length of Paired Reads\n";
	print "ARGV[3]: Error rate\n";
	print "ARGV[4]: Insert size\n";
	print "ARGV[5]: Number of reads to simulate\n";
	exit;
}
my $error_rate=$ARGV[3];
my $read_length = $ARGV[2];
my $insert_size = $ARGV[4];
my $number_of_reads = $ARGV[5];
my $acc;
my @sequence;
sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
open(FASTA,$ARGV[0]) || die("cannot open file!");
while(<FASTA>){
	my $line = $_;
	chomp($line);
	if ( $line !~ /^\>/ ){
		my @spl = split(//,$line);
		for(my $i=0; $i<scalar(@spl); $i++){
			push(@sequence,$spl[$i]);
		}
	}else{
		my @spl = split(/\>/,$line);
		$acc = $spl[1];
	}
}
close(FASTA);
my $Q = -10*log10($error_rate);
$Q += 33;
my $ascii = chr($Q);
my %index = ();
open(READ1,'>', $ARGV[1] . '.simulated' . '_' . $ARGV[2] . '_read1.fq') || die("Cannot open file!");
open(READ2,'>', $ARGV[1] . '.simulated' . '_' . $ARGV[2] . '_read2.fq') || die("Cannot open file!");
my $simulated=0;
while ( $simulated < $number_of_reads ){ 
	my $random_start = int(rand(scalar(@sequence)-2*$read_length-$insert_size));
	print READ1 "@" . "$acc" . "_" . $random_start . "_" . $simulated . "_1\n";
	print READ2 "@" . "$acc" . "_" . $random_start . "_" . $simulated . "_2\n";
	for(my $i=$random_start;$i<$random_start+$read_length;$i++){
		my $base = $sequence[$i];
		my $random_number = rand();
		if ($random_number < $error_rate){
			my $random_nuc = rand();
			if ( $base eq "A" and $random_nuc <= 0.33 ){
				print READ1 "C";
			}elsif ( $base eq "A" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ1 "G";
			}elsif ( $base eq "A" and $random_nuc >= 0.66 ){
				print READ1 "T";
			}elsif ( $base eq "C" and $random_nuc <=0.33 ){
				print READ1 "A";
			}elsif ( $base eq "C" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ1 "G";
			}elsif ( $base eq "C" and $random_nuc >= 0.66 ){
				print READ1 "T";
			}elsif ( $base eq "G" and $random_nuc <= 0.33 ){
				print READ1 "A";
			}elsif ( $base eq "G" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ1 "C";
			}elsif ( $base eq "G" and $random_nuc >= 0.66 ){
				print READ1 "T";
			}elsif ( $base eq "T" and $random_nuc <= 0.33 ){
				print READ1 "A";
			}elsif ( $base eq "T" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ1 "C";
			}elsif ( $base eq "T" and $random_nuc >= 0.66 ){
				print READ1 "G";
			}else{
				print READ1 "$base";
			}
		}else{
			print READ1 "$base";
		}
	}
	print READ1 "\n+\n";
	for (my $i=0; $i<$read_length; $i++){
		print READ1 "$ascii";
	}
	print READ1 "\n";
	for(my $i=$random_start+$read_length+$insert_size;$i<$random_start+2*$read_length+$insert_size;$i++){
		my $base = $sequence[$i];
		my $random_number = rand();
		if ($random_number < $error_rate){
			my $random_nuc = rand();
			if ( $base eq "A" and $random_nuc <= 0.33 ){
				print READ2 "C";
			}elsif ( $base eq "A" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ2 "G";
			}elsif ( $base eq "A" and $random_nuc >= 0.66 ){
				print READ2 "T";
			}elsif ( $base eq "C" and $random_nuc <=0.33 ){
				print READ2 "A";
			}elsif ( $base eq "C" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ2 "G";
			}elsif ( $base eq "C" and $random_nuc >= 0.66 ){
				print READ2 "T";
			}elsif ( $base eq "G" and $random_nuc <= 0.33 ){
				print READ2 "A";
			}elsif ( $base eq "G" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ2 "C";
			}elsif ( $base eq "G" and $random_nuc >= 0.66 ){
				print READ2 "T";
			}elsif ( $base eq "T" and $random_nuc <= 0.33 ){
				print READ2 "A";
			}elsif ( $base eq "T" and $random_nuc > 0.33 and $random_nuc < 0.66 ){
				print READ2 "C";
			}elsif ( $base eq "T" and $random_nuc >= 0.66 ){
				print READ2 "G";
			}else{
				print READ2 "$base";
			}
		}else{
			print READ2 "$base";
		}
	}
	print READ2 "\n+\n";
	for (my $i=0; $i<$read_length; $i++){
		print READ2 "$ascii";
	}
	print READ2 "\n";
	$simulated++;
}
close(READ1);
close(READ2);
