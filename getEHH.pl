#!/usr/bin/perl -w

=head1 SYNOPSIS

perl getEHH.pl [options]

 Requried parameters:
   -pairs           the space separated file containing pairs of variant
                    ids to analyse on each line with the variant with
                    the smaller genomic coordinate in the leg file
                    specified first
   -leg             the shapeit2 format legend file
   -hap             the shapeit2 format haplotype file
   -len             the length of the chromosome being analysed in bp
   -out             a label to add to the output file name
   
 Options:
   -help            brief help message
   -man             full documentation

=head1 DESCRIPTION

This program will read the given input file(s) and calculate EHH around pairs of alleles for calculating mm-iHS

 Output columns:
    Distance, in variants, from core haplotype
    
    Frequency of 00 allele combination at core haplotype
    Frequency of 01 allele combination at core haplotype
    Frequency of 10 allele combination at core haplotype
    Frequency of 11 allele combination at core haplotype
    
    EHH at position for 00 core haplotype
    EHH at position for 01 core haplotype
    EHH at position for 10 core haplotype
    EHH at position for 11 core haplotype
    
    EHH across all core haplotypes except 00
    EHH across all core haplotypes except 01
    EHH across all core haplotypes except 10
    EHH across all core haplotypes except 11
    
=cut
    
use strict;
use Tie::File;
use Getopt::Long;
use Pod::Usage;

my $man = 0;
my $help = 0;

my $input = "";
my $legfile = "";
my $hapfile = "";
my $chrlength = "";
my $label = "";
pod2usage("$0: No arguments specified.")  if ((@ARGV == 0) && (-t STDIN));
GetOptions ('pairs=s' => \$input,
            'leg=s' => \$legfile,
            'hap=s' => \$hapfile,
            'len=i' => \$chrlength,
            'out=s' => \$label,
            'help|?' => \$help, man => \$man) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;



#this file contains the list of variants that I care about. In this case it is a list of variants that cause disease and I want to see
#if particular combinations of alleles at pairs of these variants have been favoured. We do this by taking each possible combination of alleles
#at the pair of variants of interest (00, 01, 10, 11) and look at their association with the alleles at surrounding variants. For example if whenever a chromosome
#carries 00 at the pair of variants the allele at a third nearby variant is always 0 then there is a strong association between them. This script therefore
#takes each pair of variants in turn and works out how their links to surrounding alleles

print "Reading variant pairs...\n";

my @variant1 = ();
my @variant2 = ();
open(POS, $input) or die "Cannot open $input\n";
while(defined(my $line=<POS>))
{
    chomp($line);
    my @data = split(/\s+/, $line);
    push(@variant1, $data[0]);
    push(@variant2, $data[1]);
}
close POS;

#my $numPairs = scalar @pos;
#print "$numPairs pairs found...\n";


my $lineNum = 0;
my $varNum = 0;
my @snpPos = ();
my @snpId = ();
my %legArrayPos = ();
open(POS, $legfile) or die "Cannot open $legfile\n";
while(defined(my $line=<POS>))
{
    chomp $line;
    if($line !~ /^id/i)
    {
        my @data = split(/\s/, $line);
        push(@snpId, $data[0]);
        push(@snpPos, $data[1]);
        my $all = "$data[2]\t$data[3]";
        #changed this to key of id rather than pos
        $legArrayPos{$data[0]}=$lineNum;   
        $lineNum++;
    }
}
close POS;


tie my @tieArray, 'Tie::File', $hapfile or die("Unable to open file \"$hapfile\": $!\n");

#get each first SNP of interest
for(my $x = 0; $x < scalar @variant1; $x++)
{
    
    if ((!(exists($legArrayPos{$variant1[$x]}))) || (!(exists($legArrayPos{$variant2[$x]})))) {
        print "One or both variants not found in legend file  - next!\n";
        next;
    }
    my $pair1Pos = $snpPos[$legArrayPos{$variant1[$x]}];
    my $pair2Pos = $snpPos[$legArrayPos{$variant2[$x]}];
        
        
    if(abs($pair1Pos - $pair2Pos) <= 500000)
    {
        #open file unless already exists
        my $filename = "$variant1[$x]-$variant2[$x].$label.EHH.txt";
        if (-e $filename){
            print "$filename already exists - skipping...\n";
            next;
        }
        else{           
            print "Processing $filename\n";
            open(OUT, ">$filename");

            #take a window 500,000 bases either side of the SNPs of interest. Makes sure window cant extend beyond ends of the chromosomes
            my $start = $pair1Pos-500000;
            if ($start < 1) {
                $start = 1;
            }
            my $end = $pair2Pos+500000;
            if ($end > $chrlength) {
                $end = $chrlength;
            }
    
    
            #here we store the alleles for each SNP observed in each chromosome in a 2D array
            #each row in this file is a SNP and each column a chromosome (each person has two chromosomes)
            #importantly the order of the SNPs in this file is the same as in the .legend file above
            
            my $snp1ArrayPos = $legArrayPos{$variant1[$x]};
            my $snp2ArrayPos = $legArrayPos{$variant2[$x]};
            my $snp1All = "";
            my $snp2All = "";
            my @haps = ();
            my @fixed = ();
            $lineNum = 0;
            my $numHaps = 0;
            my $startArrayPos = -1;
            my $endArrayPos = -1;
            my $kill_if_SNPS_fixed = 0;
            for(my $p = 0; $p < scalar @snpPos; $p++)
            {
                #print "$start\t$end\t$snpPos[$p] - $pair1 $pair2\n";
                if ((($snpPos[$p]) >= $start) && (($snpPos[$p]) <= $end)) {
                    my $line = $tieArray[$p];
                    $haps[$lineNum] = [ split /\s/, $line ];
                    my @array = split(/\s/, $line);
                    my $allAlls = join('', @array);
                    my $c0 = $allAlls =~ tr/0//;
                    my $c1 = $allAlls =~ tr/1//;
                    if (($c0 == 0) || ($c1 == 0)) {
                        $fixed[$lineNum] =1 ;
                    }
                    else
                    {
                        $fixed[$lineNum] =0 ;
                    }
                    $numHaps = scalar @array;
                    $endArrayPos = $p;
                    if ($startArrayPos == -1) {
                        $startArrayPos = $p;
                    }
                   
                    if ($snpId[$p] eq $variant1[$x]) {
                            $kill_if_SNPS_fixed = 1 if (($c0 == 0) || ($c1 == 0)) ;
                        $snp1ArrayPos = $lineNum;
                    }
                    if ($snpId[$p] eq $variant2[$x]) {
                            $kill_if_SNPS_fixed = 1 if (($c0 == 0) || ($c1 == 0)) ;
                        $snp2ArrayPos = $lineNum;
                    }
                    $lineNum++;
                }
                last if(($snpPos[$p]) >= $end);
            }
            if ($kill_if_SNPS_fixed ==1){
                my $command = "rm $filename";
                print "One of the core variants is fixed so skipping...\n";
                `$command`;
                next;
            }

            #print header
            print OUT ("#$variant1[$x]\t$variant2[$x]\t$pair1Pos\t$pair2Pos\n");
            
            #initalise lots of arrays and hashes to store the counts
            my $num00 = 0;
            my $num01 = 0;
            my $num10 = 0;
            my $num11 = 0;
            
            my %count00 = ();
            my %count01 = ();
            my %count10 = ();
            my %count11 = ();
            
            my %countNot00 = ();
            my %countNot01 = ();
            my %countNot10 = ();
            my %countNot11 = ();
            
            my @probs00 = (0)x$snp1ArrayPos;
            my @probs01 = (0)x$snp1ArrayPos;
            my @probs10 = (0)x$snp1ArrayPos;
            my @probs11 = (0)x$snp1ArrayPos;
            
            my @probsNot00 = (0)x$snp1ArrayPos;
            my @probsNot01 = (0)x$snp1ArrayPos;
            my @probsNot10 = (0)x$snp1ArrayPos;
            my @probsNot11 = (0)x$snp1ArrayPos;
            
            #Do variants downstream of pair of gene associated variants
            #doing each chromosome one by one ($numHaps being the total number of chromosomes)
            for(my $b = 0; $b < $numHaps; $b++)
            {
                my $haplotype = "";
                #here we just count how often saw each possible combination of alleles at the pair of variants of interest
                #so if this chromsome carries a 0 at both SNPs add 1 to the 00 count
                if(($haps[$snp1ArrayPos][$b] == 0) && ($haps[$snp2ArrayPos][$b] == 0))
                {
                    $num00++;
                }
                elsif(($haps[$snp1ArrayPos][$b] == 0) && ($haps[$snp2ArrayPos][$b] == 1))
                {
                    $num01++;
                }
                elsif(($haps[$snp1ArrayPos][$b] == 1) && ($haps[$snp2ArrayPos][$b] == 0))
                {
                    $num10++;
                }
                elsif(($haps[$snp1ArrayPos][$b] == 1) && ($haps[$snp2ArrayPos][$b] == 1))
                {
                    $num11++;
                }
                
                #so we are first going down in terms of genomic cooridinates from the pair of variants of interest and seeing how haplotype identity
                #decays as distance increases
    
                #so go through each neighbouring variant in order, closest to the "left hand" variant to furthest away
                for(my $a = $snp1ArrayPos-1; $a >= 0; $a--)
                {
                    if ($fixed[$a] == 0) {
                    
                        #this stores the haplotype string. Each SNP allele is added one by one
                        #what we care about is how this string identity decays with distance. At short distances (i.e. $a is close to $snp1ArrayPos) this string
                        #should largely be the same in all chromosomes. However as distance increases recombination breaks the link between the alleles at the core SNPs
                        #and these distant SNPs so this string will increasingly become more and more different between chromosomes i.e. you will see lots of different versions
                        #of the string when the string is long
                        #so this corresponds to the following text in the voight paper:
        
                            #The EHH measures the decay of identity, as a function of distance, of haplotypes that carry a specified allele at one end. 
                            #For each allele, haplotype homozygosity [identity] starts at 1, and decays to 0 with increasing distance from the core site
        
                        #In our case we do not have a core allele though but a core combination of alleles e.g. 01 instead of a 0, as we are looking at a core pair
                        #of SNPs and not a single SNP
                        $haplotype = "$haplotype"."$haps[$a][$b]";
                        
                        #if the alleles at the two core SNPs of interest are both 0
                        if(($haps[$snp1ArrayPos][$b] == 0) && ($haps[$snp2ArrayPos][$b] == 0))
                        {
                            #here we store how often we see each version of the $haplotype string of different lengths.
                                $count00{$haplotype}++;
                           
                                $countNot01{$haplotype}++;
                                $countNot10{$haplotype}++;
                                $countNot11{$haplotype}++;
                            
                        }
                        
                        #same again but for the different allele combinations
                        elsif(($haps[$snp1ArrayPos][$b] == 0) && ($haps[$snp2ArrayPos][$b] == 1))
                        {
                                $count01{$haplotype}++;
                            
                                $countNot00{$haplotype}++;
                                $countNot10{$haplotype}++;
                                $countNot11{$haplotype}++;
                            
                        }
                        elsif(($haps[$snp1ArrayPos][$b] == 1) && ($haps[$snp2ArrayPos][$b] == 0))
                        {
                                $count10{$haplotype}++;
                            
                                $countNot01{$haplotype}++;
                                $countNot00{$haplotype}++;
                                $countNot11{$haplotype}++;
                            
                        }
                        elsif(($haps[$snp1ArrayPos][$b] == 1) && ($haps[$snp2ArrayPos][$b] == 1))
                        {
                                $count11{$haplotype}++;
                            
                                $countNot01{$haplotype}++;
                                $countNot10{$haplotype}++;
                                $countNot00{$haplotype}++;
                           
                        }
                    }
                }
            }
            
            #This is where we calculate our version of the EHH referred to in the voight paper
            #so we work out for a given distance from the core pair of SNPs, i.e. for a given string length (length($k)), what would be
            #the probability of pulling out two strings that were identical given how often we observed each string version of this length.
            #so if we only saw one version of the string at a given distance the probability would be 1.
            #Do this for each possible combination of alleles at the core pair of SNPs
        
            while ( (my $k,my $v) = each %count00 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num00)*($v/$num00);
                $probs00[$dist] = $probs00[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %count01 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num01)*($v/$num01);
                $probs01[$dist] = $probs01[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %count10 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num10)*($v/$num10);
                $probs10[$dist] = $probs10[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %count11 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num11)*($v/$num11);
                $probs11[$dist] = $probs11[$dist]+$thisProb;
            }
            
            
            #then do the same but restrict the analysis to chromosomes NOT carrying the allele combination of interest at the core pair of SNPs (analogous to iHHd in the voight paper)
        
            while ( (my $k,my $v) = each %countNot00 ) {
                my $dist = length($k)-1;
                my $numNot00 = $num01+$num10+$num11;
                my $thisProb = ($v/$numNot00)*($v/$numNot00);
                $probsNot00[$dist] = $probsNot00[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %countNot01 ) {
                my $dist = length($k)-1;
                my $numNot01 = $num00+$num10+$num11;
                my $thisProb = ($v/$numNot01)*($v/$numNot01);
                $probsNot01[$dist] = $probsNot01[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %countNot10 ) {
                my $dist = length($k)-1;
                my $numNot10 = $num00+$num01+$num11;
                my $thisProb = ($v/$numNot10)*($v/$numNot10);
                $probsNot10[$dist] = $probsNot10[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %countNot11 ) {
                my $dist = length($k)-1;
                my $numNot11 = $num00+$num01+$num10;
                my $thisProb = ($v/$numNot11)*($v/$numNot11);
                $probsNot11[$dist] = $probsNot11[$dist]+$thisProb;
            }
            
            
            #print results. The values are analogous to iHHa and iHHd values discussed in voight et al.

            for(my $a = $snp1ArrayPos-1; $a >= 0; $a--)
            {
                if ($fixed[$a] == 0) {
                    my $thisLoc = ($a*-1)-1;
                    my $arrayPos = $a;
                    print OUT ("$thisLoc\t\t$num00\t$num01\t$num10\t$num11\t\t$probs00[$arrayPos]\t$probs01[$arrayPos]\t$probs10[$arrayPos]\t$probs11[$arrayPos]\t\t$probsNot00[$arrayPos]\t$probsNot01[$arrayPos]\t$probsNot10[$arrayPos]\t$probsNot11[$arrayPos]\n");
                }
            }
            
            
            %count00 = ();
            %count01 = ();
            %count10 = ();
            %count11 = ();
            
            %countNot00 = ();
            %countNot01 = ();
            %countNot10 = ();
            %countNot11 = ();
            
            my $numValues = $lineNum-$snp2ArrayPos;
            
            @probs00 = (0)x$numValues;
            @probs01 = (0)x$numValues;
            @probs10 = (0)x$numValues;
            @probs11 = (0)x$numValues;
            
            @probsNot00 = (0)x$numValues;
            @probsNot01 = (0)x$numValues;
            @probsNot10 = (0)x$numValues;
            @probsNot11 = (0)x$numValues;
            
            #do the same but going downstream, i.e. to the "right" of the pair of variants of interest
            for(my $b = 0; $b < $numHaps; $b++)
            {
                my $haplotype = "";
                
                for(my $a = $snp2ArrayPos+1; $a < $lineNum; $a++)
                {
                    if ($fixed[$a] == 0) {
                        $haplotype = "$haplotype"."$haps[$a][$b]";
                        if(($haps[$snp1ArrayPos][$b] == 0) && ($haps[$snp2ArrayPos][$b] == 0))
                        {
                            #get the counts if match risk variant alleles
                                $count00{$haplotype}++;
                            
                                $countNot01{$haplotype}++;
                                $countNot10{$haplotype}++;
                                $countNot11{$haplotype}++;
                            
                        }
                        elsif(($haps[$snp1ArrayPos][$b] == 0) && ($haps[$snp2ArrayPos][$b] == 1))
                        {
                                $count01{$haplotype}++;
                            
                                $countNot00{$haplotype}++;
                                $countNot10{$haplotype}++;
                                $countNot11{$haplotype}++;
                            
                        }
                        elsif(($haps[$snp1ArrayPos][$b] == 1) && ($haps[$snp2ArrayPos][$b] == 0))
                        {
                                $count10{$haplotype}++;
                            
                                $countNot01{$haplotype}++;
                                $countNot00{$haplotype}++;
                                $countNot11{$haplotype}++;
                            
                        }
                        elsif(($haps[$snp1ArrayPos][$b] == 1) && ($haps[$snp2ArrayPos][$b] == 1))
                        {
                                $count11{$haplotype}++;
                            
                                $countNot01{$haplotype}++;
                                $countNot10{$haplotype}++;
                                $countNot00{$haplotype}++;
                           
                        }
                    }
                }
            }
            
            #get EHH of haplotype in question
            while ( (my $k,my $v) = each %count00 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num00)*($v/$num00);
                $probs00[$dist] = $probs00[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %count01 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num01)*($v/$num01);
                $probs01[$dist] = $probs01[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %count10 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num10)*($v/$num10);
                $probs10[$dist] = $probs10[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %count11 ) {
                my $dist = length($k)-1;
                my $thisProb = ($v/$num11)*($v/$num11);
                $probs11[$dist] = $probs11[$dist]+$thisProb;
            }
            
            #get combined EHH of others
            while ( (my $k,my $v) = each %countNot00 ) {
                my $dist = length($k)-1;
                my $numNot00 = $num01+$num10+$num11;
                my $thisProb = ($v/$numNot00)*($v/$numNot00);
                $probsNot00[$dist] = $probsNot00[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %countNot01 ) {
                my $dist = length($k)-1;
                my $numNot01 = $num00+$num10+$num11;
                my $thisProb = ($v/$numNot01)*($v/$numNot01);
                $probsNot01[$dist] = $probsNot01[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %countNot10 ) {
                my $dist = length($k)-1;
                my $numNot10 = $num00+$num01+$num11;
                my $thisProb = ($v/$numNot10)*($v/$numNot10);
                $probsNot10[$dist] = $probsNot10[$dist]+$thisProb;
            }
            while ( (my $k,my $v) = each %countNot11 ) {
                my $dist = length($k)-1;
                my $numNot11 = $num00+$num01+$num10;
                my $thisProb = ($v/$numNot11)*($v/$numNot11);
                $probsNot11[$dist] = $probsNot11[$dist]+$thisProb;
            }
            
            for(my $a = $snp2ArrayPos+1; $a < $lineNum; $a++)
            {
                my $thisLoc = $a-($snp2ArrayPos);
                my $arrayPos = $a-($snp2ArrayPos+1);
                if ($fixed[$a] == 0) {
                    print OUT ("$thisLoc\t\t$num00\t$num01\t$num10\t$num11\t\t$probs00[$arrayPos]\t$probs01[$arrayPos]\t$probs10[$arrayPos]\t$probs11[$arrayPos]\t\t$probsNot00[$arrayPos]\t$probsNot01[$arrayPos]\t$probsNot10[$arrayPos]\t$probsNot11[$arrayPos]\n");
                }                   
            }        
            close OUT;
        }
    }
}


    
