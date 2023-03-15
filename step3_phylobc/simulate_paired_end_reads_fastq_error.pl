#!/usr/bin/perl
use strict;
use warnings;
if ( scalar(@ARGV) < 6){
	print "./simulate_paired_end_reads_fastq_error.pl\n";
	print "ARGV[0]: Fasta to simulate reads from\n";
	print "ARGV[1]: Error rates per site\n";
	print "ARGV[2]: Prefix\n";
	print "ARGV[3]: Length of Paired Reads\n";
	print "ARGV[4]: Insert size\n";
	print "ARGV[5]: Number of reads to simulate\n";
	exit;
}
my $error_rate_file=$ARGV[1];
my @error_rates_per_site;
open(ERROR,$error_rate_file) || die("cannot open file!");
while(<ERROR>){
	my $line = $_;
	chomp($line);
	push(@error_rates_per_site,$line);
}
close(ERROR);
my $read_length = $ARGV[3];
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
if ( scalar(@sequence) != scalar(@error_rates_per_site) ){
	my $seqlength = scalar(@sequence);
	my $errorlength = scalar(@error_rates_per_site);
	print "Sequence length $seqlength not equal to error site length $errorlength\n";
	die;
}
my %index = ();
open(READ1,'>', $ARGV[1] . '.simulated' . '_' . $ARGV[2] . '_read1.fq') || die("Cannot open file!");
open(READ2,'>', $ARGV[1] . '.simulated' . '_' . $ARGV[2] . '_read2.fq') || die("Cannot open file!");
my $simulated=0;
while ( $simulated < $number_of_reads ){ 
	my $random_start = int(rand(scalar(@sequence)-2*$read_length-$insert_size));
	print READ1 "@" . "$acc" . "_" . $random_start . "_" . $simulated . "_1\n";
	print READ2 "@" . "$acc" . "_" . $random_start . "_" . $simulated . "_2\n";
	my @quality_scores_1;
	for(my $i=$random_start;$i<$random_start+$read_length;$i++){
		my $base = $sequence[$i];
		my $random_number = rand();
		my $Q = -10*log10($error_rates_per_site[$i]);
		$Q += 33;
		if ( $Q > 126 ){ $Q = 126; }
		my $ascii = chr($Q);
		push(@quality_scores_1,$ascii);
		if ($random_number < $error_rates_per_site[$i]){
			my $random_nuc = rand();
			if ( $base eq "A" and $random_nuc <= 0.4918 ){
				print READ1 "C";
			}elsif ( $base eq "A" and $random_nuc > 0.4918 and $random_nuc < 0.8295 ){
				print READ1 "G";
			}elsif ( $base eq "A" and $random_nuc >= 0.8295 ){
				print READ1 "T";
			}elsif ( $base eq "C" and $random_nuc <=0.5238 ){
				print READ1 "A";
			}elsif ( $base eq "C" and $random_nuc > 0.5238 and $random_nuc < 0.7899 ){
				print READ1 "G";
			}elsif ( $base eq "C" and $random_nuc >= 0.7899 ){
				print READ1 "T";
			}elsif ( $base eq "G" and $random_nuc <= 0.3754 ){
				print READ1 "A";
			}elsif ( $base eq "G" and $random_nuc > 0.3754 and $random_nuc < 0.6109 ){
				print READ1 "C";
			}elsif ( $base eq "G" and $random_nuc >= 0.6109 ){
				print READ1 "T";
			}elsif ( $base eq "T" and $random_nuc <= 0.2505 ){
				print READ1 "A";
			}elsif ( $base eq "T" and $random_nuc > 0.2505 and $random_nuc < 0.5057 ){
				print READ1 "C";
			}elsif ( $base eq "T" and $random_nuc >= 0.5057 ){
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
		print READ1 "$quality_scores_1[$i]";
	}
	print READ1 "\n";
	my @quality_scores_2;
	for(my $i=$random_start+$read_length+$insert_size;$i<$random_start+2*$read_length+$insert_size;$i++){
		my $base = $sequence[$i];
		my $random_number = rand();
		my $Q = -10*log10($error_rates_per_site[$i]);
		$Q += 33;
		if ( $Q > 126 ){ $Q = 126; }
		my $ascii = chr($Q);
		push(@quality_scores_2,$ascii);
		if ($random_number < $error_rates_per_site[$i]){
			my $random_nuc = rand();
			if ( $base eq "A" and $random_nuc <= 0.4918 ){
				print READ2 "C";
			}elsif ( $base eq "A" and $random_nuc > 0.4918 and $random_nuc < 0.8295 ){
				print READ2 "G";
			}elsif ( $base eq "A" and $random_nuc >= 0.8295 ){
				print READ2 "T";
			}elsif ( $base eq "C" and $random_nuc <=0.5238 ){
				print READ2 "A";
			}elsif ( $base eq "C" and $random_nuc > 0.5238 and $random_nuc < 0.7899 ){
				print READ2 "G";
			}elsif ( $base eq "C" and $random_nuc >= 0.7899 ){
				print READ2 "T";
			}elsif ( $base eq "G" and $random_nuc <= 0.3754 ){
				print READ2 "A";
			}elsif ( $base eq "G" and $random_nuc > 0.3754 and $random_nuc < 0.6109 ){
				print READ2 "C";
			}elsif ( $base eq "G" and $random_nuc >= 0.6109 ){
				print READ2 "T";
			}elsif ( $base eq "T" and $random_nuc <= 0.2505 ){
				print READ2 "A";
			}elsif ( $base eq "T" and $random_nuc > 0.2505 and $random_nuc < 0.5057 ){
				print READ2 "C";
			}elsif ( $base eq "T" and $random_nuc >= 0.5057 ){
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
		print READ2 "$quality_scores_2[$i]";
	}
	print READ2 "\n";
	$simulated++;
}
close(READ1);
close(READ2);
