#!/usr/bin/perl

# This program takes a simple snp file and finds the model type of  snp in windows to reduce the complexity for plotting

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTFILE, "> $ARGV[1]");

$numstrains=$ARGV[2];

# Number of columns before the genotypes

$cols=3;
# We will print out an output file that has just strain as row. 
# Columns will be sizes of LOH events
# So the length of each row will vary

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
				$length = $position[$rightside] - $position[$leftside];
				if ($length == 0)
				{
					print OUTPUTFILE "Single base recomb $chrom[$leftside] 1 leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
					$leftside=undef;
					$rightside=undef;
				}
				else 
				{
					print OUTPUTFILE "Normal $chrom[$leftside] $length leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
					$leftside=undef;
					$rightside=undef;
				}
			}
# Else basically do nothing as we are still in a non-LOH event
		}
# if B but no leftside, we do nothing
		elsif (($dataset[$b][$a] =~ /A/) | ($dataset[$b][$a] =~ /C/))
		{
# What to do if the leftside is defined (ie, we have already started an LOH), finding an A or C defines a new location
			if ($leftside) 
			{
				if ($chrom[$rightside] ne $chrom[$b])
				{
					$length = $position[$rightside] - $position[$leftside];
					if ($length == 0) {$length = 1;}
					print OUTPUTFILE "Out to telomere $chrom[$leftside] $length leftside=$leftside position=$position[$leftside] rightside=$rightside position=$position[$rightside]\t";
					$leftside=undef;
					$rightside=undef;
				}
				else
				{
					$rightside = $b;
				}
			}
			else 
			{
				$leftside = $b;
				$rightside = $b;
			}
		}
	}
	print OUTPUTFILE "\n";
}
				
						
