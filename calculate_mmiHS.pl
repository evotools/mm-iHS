#!/usr/bin/perl -w

=head1 SYNOPSIS

perl calculate_mmiHS.pl [options]

 Requried parameters:
   -files           list (with relative paths) of output files from getEHH.pl

 Options:
   -bin             the frequency bins to use when normalising (default 0.02 i.e. 2%)
   -plat            the EHH decay cutoff. Default of 0.05
   -out             a label to add to the output file name

   -help            brief help message
   -man             full documentation

=head1 DESCRIPTION

This program will read a list of files of EHH scores from getEHH.pl and calculate mm-iHS

 Output file 1 columns (.iHS.txt):
    EHH file name
    Variant 1 position
    Variant 2 position
    Alllele combination (0 = 00, 1 = 01, 2 = 10, 3 = 11)
    Standardised mm-iHS
    Frequency of core haplotype

 Output file 2 columns (.freqMeanSTDev.txt):
    Frequency bin
    Mean of unstandardised mm-iHS in bin
    Standard deviation of unstandardised mm-iHS

 Output file 3 columns (.iHH.txt):
    EHH file name
    Variant 1 position
    Variant 2 position
    Alllele combination (0 = 00, 1 = 01, 2 = 10, 3 = 11)
    iHH of allele combination
    iHH of NOT allele combination

=cut

use strict;
use POSIX qw(ceil);
use Getopt::Long;
use Pod::Usage;

my $man = 0;
my $help = 0;

my $bins = 0.02;
my $plateau = 0.05;
my $out = "out";
my $files = "";

pod2usage("$0: No arguments specified.")  if ((@ARGV == 0) && (-t STDIN));
GetOptions ('files=s' => \$files,
	    'bin=f' => \$bins,
	    'plat=f' => \$plateau,
	    'out=s' => \$out,
	    'help|?' => \$help, man => \$man) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

open(OUT, ">$out.iHS.txt") or die "Cannot open file $!";
open(OUT2, ">$out.freqMeanSTDev.txt") or die "Cannot open file $!";
open(OUT3, ">$out.iHH.txt") or die "Cannot open file $!";

print OUT "FILE\tPOS1\tPOS2\tALLELE_COMB\tSTAND_IHS\tALLELE_COMB_FREQ\n";
print OUT2 "FREQ\tMEAN\tSTD\n";

#get names of output files from previous scripts
my %files = ();
open(FILE, "$files") or die "Cannot open file $files\n $!";
while(defined(my $line2=<FILE>))
{
    chomp($line2);
    $files{$line2} = 1;
}
close FILE;


my %pos1 = ();
my %pos2 = ();
my $fileNum = 0;
my %lowest = ();

print "\n\nOpening files and filtering out truncated files and strange allele frequencies\nNumber of files checked...\n";
#go through each file i.e. each SNP pair output from previous script
while ( (my $k,my $v) = each %files )
{
    my $lineNum = 1;
    if ($fileNum % 100 == 0) {
	print "\t$fileNum\n";
    }
    $fileNum++;

    #these just keep a record of the lowest EHH value observed for each allele combination
    $lowest{$k}[0] = 1;
    $lowest{$k}[1] = 1;
    $lowest{$k}[2] = 1;
    $lowest{$k}[3] = 1;

    my $count1 = 0;
    my $count2 = 0;
    my $count3 = 0;
    my $count4 = 0;

    #open the file
    open(FILE2, "$k") or die "Cannot open file $k\n $!";
    while(defined(my $line=<FILE2>))
    {
	chomp($line);

	my @data = split(/\t+/, $line);
	if ($line !~ /^#/) {

	    #just check expected number of columns
	    if (scalar @data == 13) {

		#if value lower than previously observed lowest value for allele pair record new lowest value
		$lowest{$k}[0] = $data[5] if ($data[5] < $lowest{$k}[0]);
		$lowest{$k}[1] = $data[6] if ($data[5] < $lowest{$k}[1]);
		$lowest{$k}[2] = $data[7] if ($data[5] < $lowest{$k}[2]);
		$lowest{$k}[3] = $data[8] if ($data[5] < $lowest{$k}[3]);

		#just keep a record of how often each combination of alleles was observed i.e. so we can work out their frequency and compare to others of similar frequencies
		if($lineNum == 2) {
		    $count1 = $data[1];
		    $count2 = $data[2];
		    $count3 = $data[3];
		    $count4 = $data[4];
		}
	    }
	    else
	    {
		print "Problem - truncated file $k skipping...\n ";
		$files{$k} = 0;
	    }
	}
	#else the header line and just store SNP positions
	else
	{
	    $pos1{$k} = $data[2];
	    $pos2{$k} = $data[3];
	}
	$lineNum++;
    }
    close FILE2;
}

print "Finished filtering. Now calculating unstandardised iHS for different frequency bins.\n";
my %unStand_iHS_byFreq = ();
my %unStand_iHS_bySite = ();
my %freqs_bySite = ();
my %freq = ();
my %riskAlleles = ();
my %traits = ();
my %allf=();
my @validRegionHaps = ();
#go through each file again
while ( (my $k,my $v) = each %files )
{
    #if file read ok above
    if ($v == 1) {
	my @iHH_a = (0,0,0,0);
	my @iHH_d = (0,0,0,0);
	my @freqs = (0,0,0,0);
	my $pos1 = "";
	my $pos2 = "";
	#read file
	my $lineNum = 0;
	open(FILE2, "$k") or die "Cannot open file $k\n $!";
	while(defined(my $line=<FILE2>))
	{
	    chomp($line);
	    my @data = split(/\t+/, $line);
	    #if header
	    if ($lineNum > 0) {
		if ($lineNum  == 1) {
		    #find frequency bins
		    my $total = $data[1]+$data[2]+$data[3]+$data[4];

		    $freqs[0] = ceil(($data[1]/$total)/$bins)*$bins;
		    $freqs[1] = ceil(($data[2]/$total)/$bins)*$bins;
		    $freqs[2] = ceil(($data[3]/$total)/$bins)*$bins;
		    $freqs[3] = ceil(($data[4]/$total)/$bins)*$bins;

		    $allf{$k}{"0"} = $freqs[0];
		    $allf{$k}{"1"} = $freqs[1];
		    $allf{$k}{"2"} = $freqs[2];
		    $allf{$k}{"3"} = $freqs[3];

		    $freq{$freqs[0]} = 1;
		    $freq{$freqs[1]} = 1;
		    $freq{$freqs[2]} = 1;
		    $freq{$freqs[3]} = 1;
		}
		#for each of the four combinations of alleles
		for(my $x = 5; $x <= 8; $x++)
		{
		    #if EHH is above lowest observed value (i.e. hasnt yet plateaued) and above plateau threshold
		    if (($data[$x] > $lowest{$k}[$x-5]) && ($data[$x] > $plateau)) {
			#calculate iHH across segregating sites
			$iHH_a[$x-5] = $iHH_a[$x-5]+$data[$x];
			$iHH_d[$x-5] = $iHH_d[$x-5]+$data[$x+4];
		    }
		}
	    }
	    $lineNum++;
	}
	close FILE2;

	#this bit stores the observed unstandardised iHS in the appropriate frequency bin
	for(my $x = 0; $x < scalar @freqs; $x++)
	{
	    if (($freqs[$x] > 0) && ($iHH_d[$x] > 0)) {
		$freqs_bySite{$k}[$x] = $freqs[$x];
		my $val = "$k\t$x";
		push(@validRegionHaps, $val);
		my $iHS = log($iHH_a[$x]/$iHH_d[$x]);
		$unStand_iHS_bySite{$k}[$x] = $iHS;
		push( @{ $unStand_iHS_byFreq { $freqs[$x] } }, $iHS);
		print OUT3 ("$k\t$pos1{$k}\t$pos2{$k}\t$x\t$iHH_a[$x]\t$iHH_d[$x]\n");
	    }
	}
    }
}

print "Calculating mean and standard deviations of each frequency bin\n";
#go through each frequency bin and work out the mean and standard deviation of each
my %mean_iHS_byFreq = ();
my %stdv_iHS_byFreq = ();
while ( (my $k,my $v) = each %freq )
{
    if ($k > 0) {
	$mean_iHS_byFreq{$k} = &average(\@{ $unStand_iHS_byFreq { $k } });
	$stdv_iHS_byFreq{$k} = &stdev(\@{ $unStand_iHS_byFreq { $k } });
	print OUT2 ("$k\t$mean_iHS_byFreq{$k}\t$stdv_iHS_byFreq{$k}\n");
    }
}


#calculate and print out the standardised iHS
for(my $x = 0; $x < scalar @validRegionHaps; $x++)
{
    my @data = split("\t", $validRegionHaps[$x]);
    my $thisFreq = $freqs_bySite{$data[0]}[$data[1]];

    #just in case only one value in this frequency bin check SD isnt 0
    if ($stdv_iHS_byFreq{$thisFreq} != 0) {
	my $stand_iHS = ($unStand_iHS_bySite{$data[0]}[$data[1]]-$mean_iHS_byFreq{$thisFreq})/$stdv_iHS_byFreq{$thisFreq};
	print OUT ("$data[0]\t$pos1{$data[0]}\t$pos2{$data[0]}\t$data[1]\t$stand_iHS\t$allf{$data[0]}{$data[1]}\n");
    }

}
close OUT;
close OUT2;
close OUT3;

sub average{
	my($data) = @_;
	if (not @$data) {
		print("Empty array\n");
		return "NA";
	}
	my $total = 0;
	foreach (@$data) {
		$total += $_;
	}
	my $average = $total / @$data;
	return $average;
}
sub stdev{
	my($data) = @_;
	if(@$data == 1){
		return 0;
	}
	my $average = &average($data);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}
