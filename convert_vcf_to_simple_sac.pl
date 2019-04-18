#!/usr/bin/perl

# This program takes a vcf file and makes a file that can be used for snp analysis

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTALLELES, "> $ARGV[1]");

if (@ARGV < 2) {
    die "Improper number of arguments. Usage: perl convert_vcf_to_simple_100.pl inputvcffile outfile\n";
}

# Number of individuals in the vcf file
$numstrains=78;

# Cols is the number of columns before the real data that are variable between vcf files
$cols = 9;

# mindepth is the minimum number of reads to not filter out a pool
$mindepth = 20;

# maxdepthpool is the maximum number of reads or it should be filtered out
$maxdepth = 350;

# Threshold is the max number of strains for which missing data is possible for a locus to be included
$threshold = 20;

# Set this to have a minimum quality of a genotype to be used. The GQ score of GATK
$GQmin = 99;

# Hash of cumulative lengths
$chrtotal{chrI} = 0;
$chrtotal{chrII} = $chrtotal{chrI} + 230218;
$chrtotal{chrIII} = $chrtotal{chrII} + 813184;
$chrtotal{chrIV} = $chrtotal{chrIII} + 316620;
$chrtotal{chrV} = $chrtotal{chrIV} + 1531933;
$chrtotal{chrVI} = $chrtotal{chrV} + 576874;
$chrtotal{chrVII} = $chrtotal{chrVI} + 270161;
$chrtotal{chrVIII} = $chrtotal{chrVII} + 1090940;
$chrtotal{chrIX} = $chrtotal{chrVIII} + 562643;
$chrtotal{chrX} = $chrtotal{chrIX} + 439888;
$chrtotal{chrXI} = $chrtotal{chrX} + 745751;
$chrtotal{chrXII} = $chrtotal{chrXI} + 666816;
$chrtotal{chrXIII} = $chrtotal{chrXII} + 1078177;
$chrtotal{chrXIV} = $chrtotal{chrXIII} + 924431;
$chrtotal{chrXV} = $chrtotal{chrXIV} + 784333;
$chrtotal{chrXVI} = $chrtotal{chrXV} + 1091291;
$chrtotal{chrM} = $chrtotal{chrXVI} + 948066;


print OUTPUTALLELES "Chromosome\tPosition\tCum_Position\tSNPEFF\tGene\tAB-100-2\tAB-100-3\tAB-100-4\tAB-100-6\tAB-100-7\tAC-100-1\tAC-100-2\tAC-100-3\tAC-100-4\tAC-100-5\tAC-100-6\tAN-100-1\tAN-100-3\tAN-100-4\tAN-100-5\tAN-100-7\tAW-100-1\tAW-100-2\tAW-100-3",
"\tAW-100-4\tAW-100-6\tAY-100-1\tAY-100-2\tAY-100-3\tAY-100-4\tAY-100-5\tAY-100-6\tBB-100-3\tBB-100-5\tBB-100-6\tBB-100-7\tBB-100-8\tBC-100-3\tBC-100-4\tBC-100-5\tBC-100-6\tBC-100-7\tBC-100-8\tBN-100-2\tBN-100-3\tBN-100-4",
"\tBN-100-5\tBN-100-6\tBN-100-7\tBY-100-6\tBY-100-7\tL177_R01C1\tL178_R01C2\tL179_R01C3\tL180_R01C4\tL181_R01C5\tL182_R01C6\tL183_R01C7\tL184_R01C8\tL185_R01C9\tL186_R01C10",
"\tL187_R02C1\tL188_R02C2\tL189_R03C1\tL190_R03C2\tL191_R04C1\tL192_R04C2\tL193_R05C1\tL194_R05C2\tL195_R06C1\tL196_R06C2\tL197_R07C1\tL198_R07C2\tL199_R08C1\tL200_R08C2\tSP24\tSP245",
"\tSP253\tSP264\tSP272\tSP309-vitro\tSP309-vivo\tSP330\n";

while (<INPUTFILE>)
{
	unless ($_ =~ /#/)
	{
# Gather data
		@genos = split(/\t/, $_);
# This following line removes the whitespace at the end of the line.
		$genos[$numstrains+$cols-1] =~ s/\s+$//;
		$chr = $genos[0];
		$position = $genos[1];
# Sum up all missing genotypes to determine if locus is low quality
		$missing=0;
		$cum_position = $position + $chrtotal{$chr};
		for ($a=$cols; $a<$cols+$numstrains; $a++)
		{
			@datastrain = split(/:/, $genos[$a]);
			$genotypes[$a]= $datastrain[0];
			$depths[$a] = $datastrain[2];
			$GQ[$a] = $datastrain[3];
			if (($genotypes[$a] =~ /\.\/\./) | ($depths[$a] > $maxdepth) | ($depths[$a] < $mindepth))
			{
				$missing++;
			}
		}
# Need to take care that there could be second variants, indicted by $genos[4] =~ /\,/
		unless ($missing >= $threshold)
		{
			print OUTPUTALLELES "$chr\t$position\t$cum_position";
# Next need to get the annotation from snpEff which is stored in the 7th column
			@quals= split(/\;/, $genos[7]);
# Some column 7's have 11 items and others have 15
			if ($quals[10] =~ /ANN/)
			{
				@snpEff = split(/\|/, $quals[10]);
				print OUTPUTALLELES "\t$snpEff[1]\t$snpEff[3]";
			} 
			else 
			{
				@snpEff = split(/\|/, $quals[14]); 
				print OUTPUTALLELES "\t$snpEff[1]\t$snpEff[3]";
			}
# Next print out the genotypes based on parental and variant columns 3 and 4, respectively
			$refallele = $genos[3];
			if ($genos[4] =~ /\,/)
			{
				@parentals = split(/,/, $genos[4]);
				$parallele1 = $parentals[0];
				$parallele2 = $parentals[1];
				for ($a=$cols; $a<$cols+$numstrains; $a++)
				{
					if (($depths[$a] > $mindepth) && ($GQ[$a] >= $GQmin))
					{
						if ($genotypes[$a] =~ /0\/0/)
						{
							print OUTPUTALLELES "\t$refallele$refallele";
						}
						elsif ($genotypes[$a] =~ /0\/1/)
						{
							print OUTPUTALLELES "\t$refallele$parallele1";
						}
						elsif ($genotypes[$a] =~ /0\/2/)
						{
							print OUTPUTALLELES "\t$refallele$parallele2";
						}	
						elsif ($genotypes[$a] =~ /1\/1/)
						{
							print OUTPUTALLELES "\t$parallele1$parallele1";
						}
						elsif ($genotypes[$a] =~ /1\/2/)
						{
							print OUTPUTALLELES "\t$parallele1$parallele2";
						}
						elsif ($genotypes[$a] =~ /2\/2/)
						{
							print OUTPUTALLELES "\t$parallele2$parallele2";
						}
						else
						{
							print OUTPUTALLELES "\tNA";
						}
					}
					else
					{
						print OUTPUTALLELES "\tNA";
					}
				}
				print OUTPUTALLELES "\n";							
			}
			else
			{
				$parallele = $genos[4];
				for ($a=$cols; $a<$cols+$numstrains; $a++)
				{
					if (($depths[$a] > $mindepth) && ($GQ[$a] >= $GQmin))
					{
						if ($genotypes[$a] =~ /0\/0/)
						{
							print OUTPUTALLELES "\t$refallele$refallele";
						}
						elsif ($genotypes[$a] =~ /0\/1/)
						{
							print OUTPUTALLELES "\t$refallele$parallele";
						}
						elsif ($genotypes[$a] =~ /1\/1/)
						{
							print OUTPUTALLELES "\t$parallele$parallele";
						}
						else
						{
							print OUTPUTALLELES "\tNA";
						}
					}
					else
					{
						print OUTPUTALLELES "\tNA";
					}
				}
				print OUTPUTALLELES "\n";	
			}
		}
	}						
}			
