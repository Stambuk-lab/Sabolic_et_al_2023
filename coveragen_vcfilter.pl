#!/usr/bin/perl
use strict ; use warnings;

# coveragen_vcfilter   			# by M'Óscar 
my $version = "coveragen_vcfilter.pl";
#last edit
my $date = "26/VIII/2020";

############################

# Use to calculate mean depth per loci and prune out loci out of a specifiC range
# For options, usage and other information check the help typing the name of the program version and "help" or "--h" or so...
# coverage_vcfilter.pl -help


################################################################################
################################   PARAMETERS   ################################
# All of them can and should be set from the command line, check the help.

#my $inputname = "populations.snps.vcf";  		# input file name, should be a string either alphanumeric or alphabetic.
my $inputname = "no default";  		# input file name, should be a string either alphanumeric or alphabetic.

my $infoloci = 9;  		# check how many columns of information has the VCF file before the first sample column
#my $infocols = 9;  		# check how many columns of information has the VCF file before the first sample column

#minimum mean depth to accept
my $minwant = "min";
#my $minwant = "-4815162342";
#my $min = 4;

#maximum mean depth to accept
my $maxwant = "max";
#my $maxwant = "4815162342";
#my $max = 20;

#output file name
my $outfile = "generate_new_outname";  		#if this is not changed from here or from command line a name with all the details (populations, number of samples, number of loci), will be generated.
#my $outfile = "no_default_outname";  		#if this is not changed from here or from command line a name with all the details (populations, number of samples, number of loci), will be generated.

# do you want to prune loci?
my $noprune = 0;
#my $noprune = 0;

# ignore missing?
my $missingmatter = 0;
#my $missingmatter = 1;

#output table?
my $notable = 0;
#my $notable = 0;

#values to use if no values parsed (whole range) do no edit unless you have a stronng reason
my $lowest = -4815162342;
#my $lowest = -4815162342;
my $highest = 4815162342;
#my $highest = 4815162342;


################################################################################
################################################################################

my $helpinfo = "\n\n\t$version   $date   Help Information\n\t----------------------------------------------------\n
	This program will analyse mean coverage per locus from a vcf file.
	Genotypes need DP information encoded, and format needs to be the same for all SNPs.
	Will output a table with the mean depth per locus and can prune loci according to their mean value.
	\n\tCommand line arguments and defaults:
	--input / --vcf           Name (or path) of the VCF file. If none parsed will process all vcf files found at the working directory.
	                          If a directory is parsed will process all vcf files and save outputs there.\n
	--infocols                [int] Number of \"locus information columns\" before the first sample. Default: $infoloci\n
	--notable                 [flag] add this if you don't want the program to save a table with mean depth per locus.\n
	--noprune                 [flag] add this if you don't want the program to filter out any loci.\n
	--min                     [int] parse a minimum mean depth value per loci, any locus with mean depth below this value will be removed.\n
	--max                     [int] parse a maximum mean depth value per loci, any locus with mean depth above this value will be removed.\n
	--miss0                   [flag] add this to consider missing values as zero when calculating mean values, otherwise they will be ignored.\n\n
	Example:
	coveragen_vcfilter.pl --vcf /shared/files/datasets --noprune
	coveragen_vcfilter.pl --vcf /shared/files/datasets/trial_maf5R7maxhet6.vcf --min 2 --max 20 --notable\n\n
	This program is designed to work with the VCF files generated by the program \"populations\" (Stacks).
	One of the \"locus information columns\" from those vcf files must be \"FORMAT \", with the information about how is the DP data coded.
	in \"populations\" FORMAT column is the last of the nine \"locus information columns\", and data is coded as GT:DP:AD:GQ:GL
	other vcf file types should work fine provided that they also have a \"FORMAT\" column with \"DP\" in it.
	We don't have any relation with Stacks. This program was not tested on animals, altouth my mice look happy while I code.\n\n\n";

#Help!
my %arguments = map { $_ => 1 } @ARGV;
if(exists($arguments{"help"}) || exists($arguments{"--help"}) || exists($arguments{"-help"}) || exists($arguments{"-h"}) || exists($arguments{"--h"}) || exists($arguments{"h"})) {
	die "$helpinfo"; }










################ PASSING ARGUMENTS

use Getopt::Long;

GetOptions( "input=s" => \$inputname,    #   --input
            "vcf=s" => \$inputname,      #   --vcf
            "infocols=i" => \$infoloci,      #   --infocols
            "notable" => \$notable,      #   --notable
            "noprune" => \$noprune,      #   --noprune
            "min=i" => \$minwant,      #   --min
            "max=i" => \$maxwant,      #   --max
            "miss0" => \$missingmatter );   #   --miss0

my $welcome = "If you see this, something went wrong";

# warnings
if ($noprune == 1) {
	if ($notable == 1) { die "\n\n\tERROR!\n\n\tSo, from your command line arguments:\n\n\t--notable <- you do not want a table with the mean depth per locus\n\t--noprune <- neither want to filter-out loci according to their depth\n\n\tThen... WHY DID YOU CALL THE PROGRAM??\n\n\n\tCheck the help information, it may actually help you...\n\n$helpinfo"}
	elsif ($notable == 0) { $welcome = "calculate mean depth per locus and print it in a table"; }
}
elsif ($noprune == 0) { 
	if ($notable == 1) { $welcome = "prune loci according to their mean depth"; }
	if ($notable == 0) { $welcome = "calculate mean depth per locus, print it in a table, and prune loci according to it"; }
}

print "$version is running, it will $welcome\n";

#save values for output
my $printmax = $maxwant;
my $printmin = $minwant;
if ($maxwant eq "max" && $minwant eq "min" && $noprune == 0) { print "\n\tWARNING\n\tProgram is set to prune loci acording to their mean depth, but no accepted minimum or maximum was parsed.\n\tRidiculously low min and high max will be used to grasp all the range\n\n"; }
if ($maxwant eq "max") { $maxwant = $highest; $printmax = "max"; }
if ($minwant eq "min") { $minwant = $lowest; $printmin = "min"; }



#sort out if working with a vcf file or a path

use Cwd qw(cwd);
my $localdir = cwd;


use File::Basename;

my $singlefile = "no default";
my $workpath = "no default";
my $inputfile = "no default";

if ($inputname =~ /.*?\.vcf$/ || $inputname =~ /.*?\.VCF$/) { $singlefile = "yes"; $workpath = dirname($inputname); $inputfile = basename($inputname);}
elsif ($inputname eq "no default") { $workpath = $localdir; $singlefile = "no"; }
else { $workpath = $inputname; $singlefile = "no"; }


#define vcf filename, vcf file path to work with and directory
#check input file/s

#filter files
my @vcffiles = ();

if ($singlefile eq "no") {
	#read files
	print "No input file parsed, checking $workpath/\n";
	opendir(DIR, $workpath);						#open the directory 
	my @infiles = readdir(DIR);					#extract filenames
	closedir(DIR);
	
	foreach my $infile (@infiles) {
		next if ($infile =~ /^\.$/);				#don't use any hidden file
		next if ($infile =~ /^\.\.$/);			
		if ($infile =~ /.*?\.vcf$/) { push (@vcffiles, $infile); }		#save the vcf files
		elsif ($infile =~ /.*?\.VCF$/) { push (@vcffiles, $infile); }		#save the vcf files
	}
	
	my $filenum = scalar @vcffiles;
	if ($filenum == 0) { die "\n\n\tERROR!\n\tNo \".vcf\" or \".VCF\" files found at $workpath\nCan't work without a valid input file. Check the help information and try again... \n\n\n$helpinfo"; }
	else { $inputfile = $vcffiles[0]; } 
	
	if ($filenum > 1) { print "\n$filenum vcf files found, opening the first one...\n\n"; }
	else { $singlefile = "yes"; print "\n"; }
}
elsif ($singlefile eq "yes") { @vcffiles = ("$inputfile"); print "\n"; }






#########
#open files and calculate mean depth per locus


foreach my $thisfile (@vcffiles) {
	
	print "Reading $thisfile to calculate mean depth per locus...\n";
	#attach path in case it is in a different dir
	my $filepath = "$workpath/$thisfile"; 
	$filepath =~ s/\/\//\//g;
	$filepath =~ s/\/\//\//g;

	# open vcf file

	open my $VCFFILE, '<', $filepath or die "\nUnable to find or open $filepath: $!\n";
	
	#create output
	my %fileinfo = ();
	my %tableinfo = ('0' => "Chrom\tPos\tID\tmeanDEPTH\tsamples");
	
	
	my $k = 0;
	my $locini = 0;
	my $readformat = 8;
	my $sampcount = 0;
	my $locinum = 0;
	my $kept=0;
	my $allsamples = 0;
	
	while (<$VCFFILE>) {
		chomp;	#clean "end of line" symbols
		next if /^$/;  		#skip if blank
		next if /^\s*$/;  		#skip if only empty spaces
		my $line = $_;  		#save line
		$line =~ s/\s+$//;  		#clean white tails in lines
		my @wholeline= split('\t', $line);  		#split columns as different elements of an array
		my $numcolumn = scalar @wholeline;
		$allsamples = $numcolumn - $infoloci;
		my $lastsample = $numcolumn -1;
		my $firstsample = $infoloci;
		my $lastinfo = $infoloci - 1;
		
		if ($wholeline[0]=~/^##.*?/) { $fileinfo{$k} = $line; $k++; } #save metadata 
		elsif ($wholeline[0]=~ /^#.*?/) {
			#save index with the column that holds the format coded
			foreach my $col (0..$lastinfo) {  if ($wholeline[$col] eq "FORMAT") { $readformat = $col; }  }
			$fileinfo{$k} = $line;	#save headers
			$k++;
		}
		else {
			#check where is DP information
			my @format= split(':', $wholeline[$readformat]); #save column with the FORMAT info
			my $pos = 0;
			my $dpos = 1;
			foreach my $sigla (@format) { 
				if ($sigla eq "DP") { $dpos = $pos; } else { $pos++; } #save position of DP info
			}
			
			#calculate mean depth per locus
			my $depth = 0;
			my $divider = 0;
			#my $sampcount = 0;
			
			#loop through all samples from the row (locus)
			foreach my $col ($firstsample..$lastsample) {
				my @sample= split(':', $wholeline[$col]);  		#take one column each time
				my $isok = scalar @sample;
				my $dp = 0;
				
				
				#determine DP value and handle missing
				if ($isok == 1 && $missingmatter == 0) { $dp = -666; }
				elsif ($isok == 1 && $missingmatter == 1) { $dp = 0; }
				elsif ($missingmatter == 1 &&  $sample[$dpos] eq ".") { $dp = 0; }
				elsif ($sample[$dpos] eq ".") { $dp = -666; }
				else { $dp = $sample[$dpos]; }
				
				if ($dp != -666) { $depth = $depth + $dp; $divider++; }
				#if ($dp != -666) { $depth = $depth + $dp; $divider++; $sampcount++; }
				#elsif ($dp == -666) { $sampcount++; }
			}
			
			my $meandepth = $depth / $divider;
			
			#save table
			if ($notable == 0) { $tableinfo{$k} = "$wholeline[0]\t$wholeline[1]\t$wholeline[2]\t$meandepth\t$divider"; }
			#save good loci
			if ($noprune == 0) {  if ($meandepth >= $minwant && $meandepth <= $maxwant) { $fileinfo {$k} = $line; $kept++; }  }
			
			$k++;
			$locinum++;
			
		}
		
	}
	close $VCFFILE;
	
	print "$allsamples samples and $locinum loci analysed.\n";
	
	if ($notable == 0) {
		#file name and path
		my $outname = "no default";
		my $outpath = "no default";

		my $cleaname = $thisfile;
		$cleaname =~ s/^(.*)\.vcf/$1/;
		$cleaname =~ s/^(.*)\.VCF/$1/;
		
		$outname = "table_" . "$cleaname" . "_meanDP.txt";
		$outpath = $workpath;
		my $savepath = "$outpath/$outname";
		$savepath =~ s/\/\//\//g;
		$savepath =~ s/\/\//\//g;
		
		#print file
		open my $TOUT, '>', $savepath or die "\nUnable to create or save \"$savepath\": $!\n";
		
		#sort and print
		foreach my $line (sort {$a <=> $b} keys %tableinfo) { print $TOUT "$tableinfo{$line}\n"; }
		close $TOUT;
		print "Table with mean depth per locus saved.\n";
	}
	
	
	if ($noprune == 0) {
		
		if ($kept > 0 ) {
			#file name and path
			my $outname = "no default";
			my $outpath = "no default";

			my $cleaname = $thisfile;
			$cleaname =~ s/^(.*)\.vcf/$1/;
			$cleaname =~ s/^(.*)\.VCF/$1/;
			
			$outname = "$cleaname" . "_meanDP$printmin" . "to$printmax.vcf";
			$outpath = $workpath;
			my $savepath = "$outpath/$outname";
			$savepath =~ s/\/\//\//g;
			$savepath =~ s/\/\//\//g;
			
			#print file
			open my $FOUT, '>', $savepath or die "\nUnable to create or save \"$savepath\": $!\n";
			
			#sort and print
			foreach my $line (sort {$a <=> $b} keys %fileinfo) { print $FOUT "$fileinfo{$line}\n"; }
			close $FOUT;
			print "Saved $kept loci with mean depth equal or above $minwant and equal or below $maxwant.\n";
		}
		else { print "$kept loci found with mean depth equal or above $minwant and equal or below $maxwant.\nCheck file format and parameters"; }
	}
	
	print "\n";
}

print "\n$version is done!\n\n";



