#!/usr/bin/perl

# This program takes a simple snp file and calculates the length of the LOH tracts
# Example usage: perl calculate_LOH_lengths_sac_100_new.pl sac.LOH1.snps.only.qualfiltered.pruned.rearranged.flattened.txt sac.309.LOH.sizes_new.txt 55

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTFILE, "> $ARGV[1]");

$numstrains=$ARGV[2];

# Number of columns before the genotypes
$cols=3;

# Print out an output file that has each strain as a row. 
# Columns will be sizes of LOH events
# So the length of each row varies

$firstline = <INPUTFILE>;

@namesline = split(/\t/, $firstline);

$namesline[$numstrains+$cols] =~ s/\s+$//;

for ($a=$cols;$a<$numstrains+$cols;$a++)
{
	 $names[$a-$cols] = $namesline[$a];
}

# Below reads in the data and counts number of snps
$snp=1;

while (<INPUTFILE>)
{
	$nextline = $_;
	chomp ($nextline);
	@genos = split(/\t/, $nextline);
# This following line removes the whitespace at the end of the line.
	$genos[$numstrains+$cols] =~ s/\s+$//;
	$chrom[$snp] = $genos[0];
	$position[$snp] = $genos[1];
	for ($a=$cols;$a<$numstrains+$cols;$a++)
	{
		$dataset[$snp][$a-$cols] = $genos[$a];
	}
	$snp++;
}

$numsnps=$snp;

# Now need to parse the data and determine snp size

# Initialize values at zero
for ($a=0;$a<$numstrains;$a++)
{	
	$leftside=undef;
	$rightside=undef;
	print OUTPUTFILE "$names[$a]\t";
	for ($b=1;$b<$numsnps;$b++)
	{
		if ($dataset[$b][$a] =~ /B/)
		{
# If we already have a $leftside, then this will end the LOH event
			if ($leftside) 
			{
#				print "chrom[$leftside]=$chrom[$leftside]\t chrom[$b]=$chrom[$b]\n";
				if ($chrom[$leftside] ne $chrom[$b])
				{
					$rightside = $b-1;
					$length = $position[$rightside] - $position[$leftside];
					if ($length == 0)
					{
						print OUTPUTFILE "Chr.split.$chrom[$leftside].end leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
						$leftside=undef;
						$rightside=undef;
					}
					else 
					{
						print OUTPUTFILE "Chr ending normal $chrom[$leftside] $length leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
						$leftside=undef;
						$rightside=undef;
					}
				}
				else
				{
					$rightside = $b;
					$length = $position[$rightside] - $position[$leftside];
					print OUTPUTFILE "Normal $chrom[$leftside] $length leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
					$leftside=undef;
					$rightside=undef;
				}
			}
		}
# if B but no leftside, we do nothing
		elsif (($dataset[$b][$a] =~ /A/) | ($dataset[$b][$a] =~ /C/))
		{
# What to do if the leftside is undefined (ie, have yet to find an LOH), finding an A or C defines a new location
			unless ($leftside) 
			{
				until ($leftside)
				{
# First case is that we are at the very beginning
					if ($b-1 == 0) {$leftside = $b;}
# Second case is if the previous snp was on a different chromosome
					elsif ($chrom[$b] ne $chrom[$b-1]) {$leftside = $b;}
# Third case is if the snp is neither of the above, but we need to get to the first B, ignoring N's
					else
					{
						$previoussnp = $b-1;
						until ($leftside) 
						{
# This is a check on whether this looping led us off the end of the chromosome or the beginning of the dataset
							if (($previoussnp==0) | ($chrom[$previoussnp] ne $chrom[$b]))
							{
								$leftside=$b;
							}
							elsif ($dataset[$previoussnp][$a] =~ /B/)
							{
								$leftside = $previoussnp;
							}
							$previoussnp--;
						}
					}
				}
			}		
# Otherwise we would be in an LOH
			else
			{
# First possibility is that it jumps to another chromosome
				if ($chrom[$leftside] ne $chrom[$b])
				{
					$rightside = $b-1;
					$length = $position[$rightside] - $position[$leftside];
					if ($length == 0)
					{
						print OUTPUTFILE "Chr $chrom[$leftside] prev ending zero $length leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
						$leftside=$b;
						$rightside=undef;
					}
					else 
					{
						print OUTPUTFILE "Chr $chrom[$leftside] prev ending normal $length leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
						$leftside=$b;
						$rightside=undef;
					}
				}
			}
		}
	}
	print OUTPUTFILE "\n";
}
				
						
