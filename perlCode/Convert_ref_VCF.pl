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
	my ($chrom, $pos, $ref, $mut, $qual) = ($fields[0], $fields[1], $fields[3], $fields[4], $fields[5]);
	#my $pos2 = $pos + 1;
	print $OUTPUT "$chrom\t$pos\t$pos\t,,$qual\t$ref\t$mut\n"; 		
}
close $INPUT;
close $OUTPUT;
