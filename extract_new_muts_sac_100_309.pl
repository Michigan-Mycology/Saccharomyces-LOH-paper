#!/usr/bin/perl

# This program takes a text file from simple conversion and extracts out the new mutations
# of high quality for SSP309 (HHL)
# Essentially we look for homozygous and identical ancestral genotypes of haploid parents and diploid ancestor
# New mutations can be strains that are heterozygous or homozygous genotypes where the parents and ancestors are homozygous and different

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTALLELES, "> $ARGV[1]");

if (@ARGV < 2) {
    die "Improper number of arguments. Usage: perl extract_hq_snps_sac_100_309.pl inputtxtfile outfile\n";
}

$numstrains=78;

# Cols is the number of columns before the real data that are variable between text files
$cols = 4;

print OUTPUTALLELES "Chromosome\tPosition\tCum_Position\tSNPEFF\tGene\tAB-100-2\tAB-100-3\tAB-100-4\tAB-100-6\tAB-100-7\tAC-100-1\tAC-100-2\tAC-100-3\tAC-100-4\tAC-100-5\tAC-100-6\tAN-100-1\tAN-100-3\tAN-100-4\tAN-100-5\tAN-100-7\tAW-100-1\tAW-100-2\tAW-100-3",
"\tAW-100-4\tAW-100-6\tAY-100-1\tAY-100-2\tAY-100-3\tAY-100-4\tAY-100-5\tAY-100-6",
"\tL177_R01C1\tL178_R01C2\tL179_R01C3\tL180_R01C4\tL181_R01C5\tL182_R01C6\tL183_R01C7\tL184_R01C8\tL185_R01C9\tL186_R01C10",
"\tL187_R02C1\tL188_R02C2\tL189_R03C1\tL190_R03C2\tL191_R04C1\tL192_R04C2\tL193_R05C1\tL194_R05C2\tL195_R06C1\tL196_R06C2\tL197_R07C1\tL198_R07C2\tL199_R08C1\tL200_R08C2\tSP24",
"\tSP253\tSP309-vitro\tSP309-vivo\tNumhets\n";

while (<INPUTFILE>)
{
	chomp($_);
	@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
	$genos[$numstrains+$cols-1] =~ s/\s+$//;
# First check to see if all of the haploids and their derived diploids are the same
# SP24 (75), SP253 (77), 309-vivo (80), 309-vitro (81) 
	if (($genos[75] eq $genos[77]) && ($genos[75] eq $genos[80]) && ($genos[75] eq $genos[81]))
	{
		print "\nfirst round";
# Make sure they are not NA or heterozyg, just check for homozygosity, one check should do for all
		if (($genos[75] =~ /AA/) | ($genos[75] =~ /CC/) | ($genos[75] =~ /GG/) | ($genos[75] =~ /TT/))
		{
# Then check each of the evolved. Look to see if they are different than 79 and not NA
			$mutation= "F";
			print "\ninside";
			for ($a=5;$a<32;$a++)
			{			
				if (($genos[$a] ne $genos[81]) && ($genos[$a] ne 'NA'))
				{
					$mutation = "T";
					print "\nsecond round";
				}
			}
			for ($a=51;$a<75;$a++)
			{			
				if (($genos[$a] ne $genos[81]) && ($genos[$a] ne 'NA'))
				{
					$mutation = "T";
				}
			}
			if ($mutation eq "T")
			{
				$hetsum=0;
				print OUTPUTALLELES "$genos[0]\t$genos[1]\t$genos[2]\t$genos[3]\t$genos[4]\t";
				for ($a=5;$a<32;$a++)
				{
					if (($genos[$a] ne 'AA') && ($genos[$a] ne 'CC') && ($genos[$a] ne 'GG') && ($genos[$a] ne 'TT') && ($genos[$a] ne 'NA'))
					{
						$hetsum++;
					}
					print OUTPUTALLELES "$genos[$a]\t";
				}
				for ($a=51;$a<75;$a++)
				{
					if (($genos[$a] ne 'AA') && ($genos[$a] ne 'CC') && ($genos[$a] ne 'GG') && ($genos[$a] ne 'TT') && ($genos[$a] ne 'NA'))
					{
						$hetsum++;
					}
					print OUTPUTALLELES "$genos[$a]\t";
				}
				print OUTPUTALLELES "$genos[75]\t$genos[77]\t$genos[80]\t$genos[81]\t$hetsum\n";
			}
		}
	}
}
		  
