#!/usr/bin/perl

# This program takes a text file from simple conversion and extracts out the new mutations
# of high quality for SSP272 (LHL)
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

print OUTPUTALLELES "Chromosome\tPosition\tCum_Position\tSNPEFF\tGene\tBB-100-3\tBB-100-5\tBB-100-6\tBB-100-7\tBB-100-8\tBC-100-3\tBC-100-4\tBC-100-5\tBC-100-6",
"\tBC-100-7\tBC-100-8\tBN-100-2\tBN-100-3\tBN-100-4\tBN-100-5\tBN-100-6\tBN-100-7\tBY-100-6\tBY-100-7\tSSP245\tSSP264\tSSP272\tNumhets\n";

while (<INPUTFILE>)
{
	chomp($_);
	@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
	$genos[$numstrains+$cols-1] =~ s/\s+$//;
# First check to see if all of the haploids and their derived diploids are the same
# SP245 (76), SP264 (78), SP272 (79) 

	if (($genos[76] eq $genos[78]) && ($genos[76] eq $genos[79]))
	{
# Make sure they are not NA or heterozyg, just check for homozygosity, one check should do for all
		if (($genos[76] =~ /AA/) | ($genos[76] =~ /CC/) | ($genos[76] =~ /GG/) | ($genos[76] =~ /TT/))
		{
# Then check each of the evolved. Look to see if they are different than 79 and not NA
			$mutation= "F";
			for ($a=32;$a<51;$a++)
			{			
				if (($genos[$a] ne $genos[79]) && ($genos[$a] ne 'NA'))
				{
					$mutation = "T";
				}
			}
			if ($mutation eq "T")
			{
				$hetsum=0;
				print "should print\n";
				print OUTPUTALLELES "$genos[0]\t$genos[1]\t$genos[2]\t$genos[3]\t$genos[4]\t";
				for ($a=32;$a<51;$a++)
				{
					if (($genos[$a] ne 'AA') && ($genos[$a] ne 'CC') && ($genos[$a] ne 'GG') && ($genos[$a] ne 'TT') && ($genos[$a] ne 'NA'))
					{
						$hetsum++;
					}
					print OUTPUTALLELES "$genos[$a]\t";
				}
				print OUTPUTALLELES "$genos[76]\t$genos[78]\t$genos[79]\t$hetsum\n";
			}
		}
	}
}
		  
