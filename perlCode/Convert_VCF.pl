use warnings;
use strict;

if (@ARGV != 2) {
	die "need to provide 2 inputs: shared variant file and outputfile name\n";
}

my ($inputfile, $outputfile) = @ARGV;

open (my $INPUT, "<", $inputfile);
open (my $OUTPUT, ">", $outputfile);

while(<$INPUT>) {
	chomp;
	my $line = $_;
	next if ($line =~ m/^\#/);
	my @fields = split;
	my ($chrom, $pos, $ref, $mut, $qual, $info) = ($fields[0], $fields[1], $fields[3], $fields[4], $fields[5], $fields[9]);
	my @splitinfo = split(/\:/, $info);
	my ($refnum, $mutnum, $varfreq);
	my @nucnums = split(/\,/, $splitinfo[1]);
	($refnum, $mutnum) = ($nucnums[0], $nucnums[1]);
	my $totalnum = $refnum + $mutnum;
	if ($totalnum == 0){
	$varfreq = "NAN";
	} else {
	$varfreq = sprintf("%.3f",$mutnum/$totalnum);
	}
	print $OUTPUT "$chrom\t$pos\t$pos\t$totalnum,$mutnum,$qual\t$ref\t$mut\t$varfreq\n"; 		
}
close $INPUT;
close $OUTPUT;
