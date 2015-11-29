#!/usr/bin/perl -w
use strict;
use Getopt::Long;


###############################################################
##################### DOCUMENTATION ###########################

## Author: Sayantan Das
## Date: July 2014
## Goal: Grep a file (with one column) in another file (will grep in the first column of target file
my $location="/net/fantasia/home/sayantan/IMPUTATION/Test_MinimacX_Take_2/test";
my $title = " [ConformRefAltVCF]";

# print "$title\n";
# print "\t\@: $location\n";
# print "\tby: Sayantan Das (sayantan\@umich.edu)\n\n";
 
## Input Parameters
## 

## Input Files
## (1) f: File to search for (with only one column)
## (2) t: File in which to search for (can of course have more than one column, will search in first column)
## (3) o: Output File Name

## Output Files

################# END OF DOCUMENTATION ########################
###############################################################

my $MIN_NUM_OPTS = 3;

if($#ARGV<($MIN_NUM_OPTS*2 - 1))
{
	&usage();
	die "Failed to parse parameters\n\n";
}


my %opts = ();
$opts{c}=1;

Getopt::Long::GetOptions(\%opts,qw(
	f=s
	t=s
	c=s
	o=s
)) || die "Failed to parse options\n\n";


&printPar();

my %dosesnps = (); # key: marker name; # value: index (0,...);
my @doseIndices = (); # indices (0,..) in .snps for SNPs in both .dose and .ped;
my @pedIndices = (); # indices (0,..) in .dat for SNPs in both .dose and .ped;
my @snps = ();
my @majAlleles = ();
#my @refAlleles = ();
my %refAlleles = ();
my %refVCFsnps = ();
my %refVCFref = ();
my %refVCFalt = ();
my %goldsnpLists = ();
my %doseIDs = (); # key: ID; #vlaue: line index (0,..) in .dose;
my %goldIDs = (); # key: ID; #vlaue: line index (0,..) in .dose;
my %relevantPedLines = (); # key: ID; #value: entire line;
my @sumxy = ();
my @sumx = ();
my @sumy = ();
my @goldNames=();
my @sumsqx = ();
my @sumsqy = ();
my @numObsG = ();
my @doseR2 = ();

my %ACGTcodes = ( A => 1, C => 2, G => 3, T => 4, 1 => 1, 2 => 2, 3 => 3, 4 => 4);

if ($opts{f} =~ /(.+)\.gz$/) { $opts{f} = $1; }   # Remove gz suffix if provided
if ($opts{t} =~ /(.+)\.gz$/) { $opts{t} = $1; }   # Remove gz suffix if provided

#if ($opts{ped} =~ /(.+)\.gz$/) { $opts{ped} = $1; }   # Remove gz suffix if provided

&main(); 


# begin sub

sub usage
{
	# print "\n";
	# print "Usage:\n\t";
	# print "-f \t File to search  (required) \n\t";
	# print "-t \t File in which to search (required) \n\t";
	# print "-c \t Column of -t to seach in (default=1) \n\t";
	# print "-o \t Output file name (Required) \n\t";
	# print "\n";
}

sub printPar
{
	# print "\n";
	# print "Parameters in Effect:\n";
	# print "\t File to search  \n\t\t-f '$opts{f}'\n";
	# print "\t File in which to search \n\t\t-t '$opts{t}'\n";
	# print "\t Column of -t to seach in \n\t\t-c '$opts{c}'\n";
	# print "\t Output file name \n\t\t-o '$opts{o}'\n";
	# print "\n";
}



sub rm
{
	if(-e $_[0])
	{
		# print "WARNING: $_[0] existed and removed!!\n\n";
		system "/bin/rm -rf $_[0]";
	}
}



sub bakrm
{
	if(-e $_[0])
	{
		my $tt = $_[0].".old";
		# print "WARNING: $_[0] existed and moved to $tt !!\n\n";
		system "mv $_[0] $tt";
	}
}


sub swap
{
	my $tmp = $_[0];
	$_[0] = $_[1];
	$_[1] = $tmp;
}


sub main
{
	&grep();
}



sub grep
{
	# (1) on markers;
	my $linenum = 0;
	my $size=0;
	my %items = ();
	my $ffile;
	my $tarfile;
	
		
	my $ffilegz = $opts{f}.".gz";
	if (-r $ffilegz)
	{
		open($ffile, "gunzip -c $ffilegz |") || die "Can't open compressed file $ffilegz\n\n";
	}
	else
	{
		open($ffile, $opts{f}) || die "Can't openc file $opts{f}\n\n";
	}

	while(my $line = <$ffile>)
	{
		$linenum++;
		chomp($line);
		$line =~ s/^\s+//;
		my @lineArr = (split(/\s+/,$line));
		if ($#lineArr>0)
		{
			print "-f file [$opts{f}] cannot have more than one column (at line $linenum) \n\n";
			die;
		}
		if (!exists $items{$lineArr[0]})
		{
			$items{$lineArr[0]}=$linenum;
		}
		else
		{
			print "-f file cannot have duplicates [$opts{f}]\n";
			print "Duplicate Entry is $lineArr[0] \n";
			die;
		}
		
	}
	# print "\n -f file [$opts{f}] has $linenum entries \n\n";
	close($ffile)	;
	
	
	my $tfilegz = $opts{t}.".gz";
	if (-r $tfilegz)
	{
		open($tarfile, "gunzip -c $tfilegz |") || die "Can't open compressed file $tfilegz\n\n";
	} 
	else
	{
		open($tarfile, $opts{t}) || die "Can't openc file $opts{t}\n\n";
	}
	
	my $of = $opts{o};
	my $miss = $opts{o}."miss";
	&rm($of);
	&rm($miss);
	
	open(OUT,">>$of") || die "Can't append file $of\n\n";
	open(MISS,">>$miss") || die "Can't append file $miss\n\n";
	
	my $count=0;
	my $nomatch=0;
	$linenum=0;
	while(my $line = <$tarfile>)
	{
		$linenum++;
		chomp($line);
		$line =~ s/^\s+//;
		my @lineArr = (split(/\s+/,$line));
		
		if (exists $items{$lineArr[$opts{c}-1]})
		{
			print OUT "$line\n";
			$count++;
		}
		else
		{
			
			print MISS "$line\n";
			
			$nomatch=$nomatch+1;
		}
		
	}
	close($tarfile);
	close($of);
	
	# print "\n $count entries extracted from target file \n\n";
	# print "\n NOWWWW sorting them ...";
	 # system("awk 'NR==FNR{o[FNR]=\$2; next} {t[\$1]=\$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' $opts{f} $opts{o} > $opts{o}.temp");
 #system("rm $opts{o}");

			
	}

	
	
	

