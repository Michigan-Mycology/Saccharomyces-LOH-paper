#!/usr/bin/perl

# This program takes a vcf file and makes a text file with quality measure of each snp
# This file can then be analyzed in R to look at the distribution of measures

# open input and output files
open (INPUTFILE, $ARGV[0]);
open (OUTPUTQC, "> $ARGV[1]");

print OUTPUTQC "Contig\tPosition\tQUAL\tAC\tAN\tBaseQRankSum\tClippingRankSum\tFS\tDP\tMQ\tMQRankSum\tQD\tReadPosRankSum\tSOR\n";

while (<INPUTFILE>)
{
	unless ($_ =~ /#/)
	{
			@genos = split(/\t/, $_);
			unless (($genos[4] =~ /\*/) | ($genos[4] =~ /\,/))
			{
				print OUTPUTQC "$genos[0]\t$genos[1]\t$genos[5]\t";
				$QCdata = $genos[7];
				my $AC = $QCdata =~ /AC=([a-z0-9-.]+)/;
				print $genos[0], " ", $genos[1], " ", $1,"\n";
				if ($AC) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}
				my $AN = $QCdata =~ /AN=([a-z0-9-.]+)/;
				if ($AN) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}
				my $BQRS = $QCdata =~ /BaseQRankSum=([a-z0-9-.]+)/;
				if ($BQRS) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}
				my $CRS = $QCdata =~ /ClippingRankSum=([a-z0-9-.]+)/;
				if ($CRS) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}				
				my $FS = $QCdata =~ /FS=([a-z0-9-.]+)/;
				if ($FS) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}				
				my $DP = $QCdata =~ /DP=([a-z0-9-.]+)/;
				if ($DP) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}				
				my $MQ = $QCdata =~ /MQ=([a-z0-9-.]+)/;
				if ($MQ) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}				
				my $MQRS = $QCdata =~ /MQRankSum=([a-z0-9-.]+)/;
				if ($MQRS) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}				
				my $QD = $QCdata =~ /QD=([a-z0-9-.]+)/;
				if ($QD) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}		
				my $RPRS = $QCdata =~ /ReadPosRankSum=([a-z0-9-.]+)/;
				if ($RPRS) {print OUTPUTQC "$1\t";}
				else {print OUTPUTQC "NA\t";}			
				my $SOR = $QCdata =~ /SOR=([a-z0-9-.]+)/;
				if ($SOR) {print OUTPUTQC "$1\n";}
				else {print OUTPUTQC "NA\n";}			
			}
	}
}
