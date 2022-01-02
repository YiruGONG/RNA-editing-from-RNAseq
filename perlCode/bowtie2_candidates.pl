#BLAT mismatched reads to ensure unique mapping!

use warnings;
use strict;

if (@ARGV != 5) {
	die "need to provide 5 input: Core number, Variant list, INDEXED BAM alignment file, softmasked reference genome, and output file name\n";
}
my ($ncore, $inputfile, $bamfile, $refGenome, $outputfile) = ($ARGV[0], $ARGV[1], $ARGV[2],$ARGV[3],$ARGV[4]);

my $minbasequal = 25;
my $minmismatch = 1;
my $scorelimit = 0.95; #second blat hit score must be less than 95% of first hit!


open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);

my $fafile = join '', $inputfile, '.fatmp';
my $samfile = join '', $inputfile, '.samtmp';
open (FAFILE, ">", $fafile);

while (<$INPUT>) {
	chomp;
	my $inputline = $_;
	my @fields = split;
	my $TEMPNAME = join '', $outputfile,'_tmp';
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $bamposition = join '', $chr,':',$position,'-',$position;
	system("samtools view $bamfile $bamposition > $TEMPNAME"); #get reads overlapping candidate site
	my $editnuc = $fields[5];
	my $newmismatch = 0;
	my $mismatchreadcount = 0;

	open(my $TMPFILE, "<", $TEMPNAME);
	while (<$TMPFILE>) {
		chomp;
		my @bamfields = split;
		my ($alignment, $readstart, $cigar, $sequence, $qualities) = ($bamfields[1], $bamfields[3], $bamfields[5], $bamfields[9], $bamfields[10]);
		my @sequencebases = split(//,$sequence);
		my @qualscores = split(//,$qualities);
		my ($currentpos, $readpos) = ($readstart,1);
		my $base_readpos;
		my @cigarnums = split(/[MIDNSHP]/, $cigar);
		my @cigarletters = split(/[0-9]+/,$cigar);
		shift @cigarletters;

		for (my $i = 0; $i < @cigarnums; $i++) {
			if ($cigarletters[$i] =~ m/[ISH]/) {
				$readpos = $readpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/[DN]/) {
				$currentpos = $currentpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/M/) {
				for (my $j = 0; $j < $cigarnums[$i]; $j++) { #for each read, determine if it contains the mismatch at the candidate site
					if ($currentpos -$readstart <  @sequencebases && $readpos < @sequencebases){
						$base_readpos = 1 if (($currentpos == $position) && ($sequencebases[$readpos-1] eq $editnuc) && (ord($qualscores[$readpos-1]) >= $minbasequal+33)); #SANGER format!
					}
					$currentpos++;
					$readpos++;	
				}	
			}
		}
		if ($base_readpos) {
			print FAFILE ">$chr\-$position\-$mismatchreadcount\n$sequence\n"; #make fasta file of all mismatched reads
			$mismatchreadcount++;
		}
	}
	close $TMPFILE;
	system("rm $TEMPNAME");
}
close $INPUT;


system("bowtie2 -p $ncore --local -a -f --rfg 0,0 --score-min G,20,0 --no-hd -x $refGenome -U $fafile -S $samfile");
#	--ma 1 --mp 0,0 --np 0 --rfg 0,0 --rdg 0,0 ##test qualscore changes (qual and seq part)
#--rfg 0,0 --rdg 0,0 : no gap penalty
system("rm $fafile");



open(SAM, "<", $samfile );
my %samhash;
while(<SAM>) { #parse SAM outputfile: for each mismatched read that was blatted, hash together all the blat hits 
	chomp;
	my @samfields = split;
	my ($name,$readstart, $cigar) = ($samfields[0],$samfields[3],$samfields[5]); #qName: query sequence name chr_pos
	my ($currentpos, $readpos, $blockStarts) = ($readstart,1, $readstart);
	my $base_readpos,$blockSizes;
	my ($matchnum, $blockCount) = (0,1);
	my @cigarnums = split(/[MIDNSHP]/, $cigar);
	my @cigarletters = split(/[0-9]+/,$cigar);
	shift @cigarletters;
	
	for (my $j = 0; $j < @cigarnums; $j++) {
		if ($cigarletters[$j] =~ m/[M]/) {
			$matchnum= $matchnum + $cigarnums[$j];
			$currentpos = $currentpos + $cigarnums[$j];
			$readpos = $readpos + $cigarnums[$j];
			#可判断是否对应
		}
		elsif ($cigarletters[$j] =~ m/[DN]/) {
			$blockSizes = join ",", $blockSizes, $cigarnums[$j]);
			$blockStarts = join ",", $blockStarts, $currentpos;
			$blockCount++;
			$currentpos = $currentpos + $cigarnums[$j];
		}
		elseif ($cigarletters[$j] =~ m/[ISH]/) {
			$readpos = $readpos + $cigarnums[$j];
		}
	}
	$blockSizes = $readpos unless ($blockSizes);

	my $btscore = join '@',$matchnum,$samfields[2],$blockCount,$blockSizes,$blockStarts; #matches, tName, blockCount, blockSizes, tStarts
	if ($samhash{$name}) {
		$samhash{$name} = join '-', $samhash{$name}, $btscore;
	} elsif ( !$samhash{$name}) {
		$samhash{$name} = $btscore; 
	}
}
close SAM;	

my %sitehash;
my %discardhash;
foreach my $samkey (keys %samhash) { #look at each blatted read one by one
	my @splitkey = split(/\-/, $samkey);
	my $site = join '_', $splitkey[0],$splitkey[1]; #assign this read to its corresponding candidate site
	my @samlines = split(/\-/, $samhash{$samkey}); #separate out the blat hits for each read
	my $largestscore = 0;
	my $largestscoreline = $samlines[0];
	my @scorearray;
	foreach my $scoreline (@samlines) { #find the largest scored blat hit in the genome
		my @scoresarray = split(/\@/,$scoreline);
		my $linescore = $scoresarray[0];
		push(@scorearray,$linescore);
		if ($linescore > $largestscore) {
			$largestscoreline = $scoreline;
			$largestscore = $linescore;
		}
	}
	@scorearray = sort {$b <=> $a} @scorearray; #sort the blat hits by score
	$scorearray[1] = 0 unless ($scorearray[1]);
	my @splitlargestline = split(/\@/,$largestscoreline);
	my $overlapfound = 0;
	if ($splitlargestline[1] eq $splitkey[0] && $scorearray[1] < ($scorearray[0]*$scorelimit)) { #ensure that score of second blat hit is less than 95% of first hit
		my ($numblocks, $blocksizes, $blockstarts) = ($splitlargestline[2],$splitlargestline[3],$splitlargestline[4]);
		my @blocksizes = split(/\,/,$blocksizes);
		my @blockstarts = split(/\,/,$blockstarts);
		for (my $i = 0; $i < $numblocks; $i++) {
			my $startpos = $blockstarts[$i]+1;
			my $endpos = $blockstarts[$i] + $blocksizes[$i];
			$overlapfound = 1 if ($splitkey[1] >= $startpos && $splitkey[1] <= $endpos); #ensure that first blat hit overlaps the candidate editing site
		}
		if ($sitehash{$site} && $overlapfound) { #check if this read passes the blat criteria
			$sitehash{$site}++;
		} elsif (!$sitehash{$site} && $overlapfound) {
			$sitehash{$site} = 1;
		}
	}
	unless ($overlapfound) { #check if this read has failed the blat criteria
		if ($discardhash{$site}) {
			$discardhash{$site}++;
		} elsif (!$discardhash{$site}) {
			$discardhash{$site} = 1;
		}
	}
}

open (SITES2, "<", $inputfile ) or die "error opening inputfile: $!\n"; #open input file again and check if each candidate passes the blat filtering
while(<SITES2>) {
	chomp;
	my @fields = split;
	my $inputline = $_;
	my ($cov,$oldalter) = split(/\,/,$fields[3]);
	my $name = join '', $fields[0],'_',$fields[1];
	if ($sitehash{$name}) {
		my $newalter = $sitehash{$name};
		my $discardnum;
		if ($discardhash{$name}) {
			$discardnum = $discardhash{$name};
		} else {
			$discardnum = 0;
		}
		my $newcov = $cov - ($oldalter - $newalter);
		my $neweditfreq = sprintf("%.3f", $newalter/$newcov);
		print $OUTPUT "$fields[0]\t$fields[1]\t$fields[2]\t$newcov,$newalter\t$fields[4]\t$fields[5]\t$neweditfreq\n" if ($newalter >= $minmismatch && $newalter > $discardnum); #make sure that number of read passing blat criteria is greater than number that fail
	}
}
close SITES2;
close $OUTPUT;
system("rm $samfile");
