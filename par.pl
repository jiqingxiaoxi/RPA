use strict; use warnings;
use Getopt::Long;

##used in read parameters
my $prefix;
my $common_file;
my $specific_file;
my $par_PrimerCommon=0;
my $par_PrimerSpecific=4;
my $par_ProbeCommon=0;
my $par_ProbeSpecific=3;
my $flag;
my $help;
my $path_bowtie2;
my $index;
my $NoProbe=0;
my $dir;
my $left;
my $check;
my $thread;
my $ref;
my $seq="";
my $value;
my $begin;
my $stop;

my $line;
my @array;
my $pos;
my $len;
my $i;
my $primer;
my $name;
my $mismatch;
my @list_index;
my %common;
my %specific;
my @list_common;
my @file;
my @status;
my $start;
my $end;
my $take;
my $num_common=0;
my $num_specific=0;
my $MisPos;
my $MD;

$help="USAGE:
  perl $0 --in <input_name> --ref <ref_genome> --common[--specific] <genomes_list> --bowtie2 <bowtie2> --index <indexs_bowtie2> [options]*\n
ARGUMENTS:
  --in <input_name>
    the file name of candidate single primer/probe, files are generated from Single program
  --ref <ref_genome>
    reference genome, fasta formate
  --dir <directory>
    dirctory for files of candidate single primer/probe regions
    default: current directory
  --NoProbe
    without exo probe
  --common <genomes_list>
    the genomes in the file(target genomes) are expected to be amplified by RPA primer sets
  --specific <genomes_list>
    the genomes in the file(background genomes) are not expected to be amplified by RPA primer sets
  --left
    background_group = all_genome_in_database - target_group
    used with --common
    invalid if exist --specific
  --bowtie2 <bowtie2>
    the bowtie2 program
  --index <database>
    bowtie2 index file name, comma-separated
  --Primer_s <int>
    the max number of mismatches allowed when align single primers to background genomes
    the bigger of the value, the more specific
    default: 4
  --Primer_c <int>
    the max number of mismatches allowed when align single primers to target genomes
    the smaller of the value, the more common
    default: 0
  --Probe_s <int>
    the max number of mismatches allowed when align probe to background genomes
    the bigger of the value, the more specific
    default: 3
  --Probe_c <int>
    the max number of mismatches allowed when align probe to target genomes
    the smaller of the value, the more common
    default: 0
  --threads <int>
    number of threads to launch when align
    default: 1
  --help|--h
    print help information\n";

$check=GetOptions("in=s"=>\$prefix,"ref=s"=>\$ref,"dir=s"=>\$dir,"common=s"=>\$common_file,"specific=s"=>\$specific_file,"left"=>\$left,"bowtie2=s"=>\$path_bowtie2,"index=s"=>\$index,"Primer_c=i"=>\$par_PrimerCommon,"Primer_s=i"=>\$par_PrimerSpecific,"Probe_c=i"=>\$par_ProbeCommon,"Probe_s=i"=>\$par_ProbeSpecific,"threads=i"=>\$thread,"help|h"=>\$flag,"NoProbe"=>\$NoProbe);

if($check==0)
{
	print $help;
	exit;
}
##check input paramters
if($flag)
{
        print $help;
        exit;
}
if(!$prefix)
{
	print "Error!\nDon't have the --in! The name of candidate primer/probe files should be supplied.\n";
	print $help;
	exit;
}

if(!$ref)
{
        print "Error!\nDon't have the --ref! The reference genome file should be supplied.\n";
        print $help;
        exit;
}

if(!$dir)
{
	$dir=$ENV{'PWD'};
	$dir=$dir."/";
}
else
{
	if($dir!~/\/$/)
	{
		$dir=$dir."/";
	}
}

if((!$common_file) && (!$specific_file) && (!$left))
{
	print "Error!\nOne of --common, --specific and --left should be existed.\n";
	exit;
}

if(!$path_bowtie2)
{
	print "Error!\nDon't have the --bowtie2! The bowtie program should be supplied.\n";
	print $help;
	exit;
}
else
{
	if(!(-e $path_bowtie2))
	{
		print "Don't have the bowtie2 program in $path_bowtie2!\n";
		exit;
	}
}

if(!$index)
{
	print "Error!\nDon't have --index! Comma-separated index files used in bowtie2.\n";
	print $help;
	exit;
}
@list_index=split(",",$index);

if($par_PrimerCommon>$par_PrimerSpecific)
{
	print "Error!\nThe argument of --Primer_c (now is $par_PrimerCommon) must <= the argument of --Primer_s (now is $par_PrimerSpecific)!\n";
	exit;
}
if(($NoProbe==0)&&($par_ProbeCommon>$par_ProbeSpecific))
{
	print "Error!\nThe argument of --Probe_c (now is $par_ProbeCommon) must <= the argument of --Probe_s (now is $par_ProbeSpecific)!\n";
	exit;
}

if(!$thread)
{
	$thread=1;
}

###################################prepare
if($common_file)
{
	open(IN,"<$common_file") or die "Can't open the $common_file inputed from --common!\n";
	while(<IN>)
	{
		chomp;
		$line=$_;
		if($line=~/^\>/)
		{
			($line)=$line=~/^\>(.+)$/;
		}
		@array=split(" ",$line);
		if(length($array[0])<=300)
		{
			if(exists $common{$array[0]})
			{
				print "Warnings: in $common_file, the name of \"$array[0]\" exists in different genomes.\n";
				next;
			}
			$common{$array[0]}=$num_common;
			$list_common[$num_common]=$array[0];
			$num_common++;
		}
		else
		{
			$flag=substr($array[0],0,300);
			if(exists $common{$flag})
			{
				print "Warnings: in $common_file, the name of \"$flag\" exists in different genomes.\n";
				next;
			}
			$common{$flag}=$num_common;
			$list_common[$num_common]=$flag;
			$num_common++;
		}
	}
	close IN;
}

if($specific_file)
{
	open(IN,"<$specific_file") or die "Can't open the $specific_file file inputed from --specific!\n";
	while(<IN>)
	{
		chomp;
		$line=$_;
		if($line=~/^\>/)
                {
                        ($line)=$line=~/^\>(.+)$/;
                }
                @array=split(" ",$line);
		if(length($array[0])<=300)
		{
			if(exists $specific{$array[0]})
			{
				print "Warnings: in $specific_file, the name of \"$array[0]\" exists in different genomes.\n";
				next;
			}
			$specific{$array[0]}=$num_specific;
			$num_specific++
		}
		else
		{
			$flag=substr($array[0],0,300);
			if(exists $specific{$flag})
			{
				print "Warnings: in $specific_file, the name of \"$flag\" exists in different genomes.\n";
				next;
			}
			$specific{$flag}=$num_specific;
			$num_specific++;
		}
	}
	close IN;
}

$flag=0;
open(IN,"<$ref") or die "Can't open $ref file!\n";
while(<IN>)
{
        chomp;
        $line=$_;
        if($line=~/^\>/)
        {
		$flag++;
                next;
        }
        $seq=$seq.$line;
}
close IN;
if($flag>1)
{
	print "Error!The $ref file contains more than one sequences.\n";
	exit;
}
#####################################run bowtie2 and analysis for primers
print "Now the program is handling the single primer file\n";
$start=time();
$file[0]=$dir."Primer/".$prefix;
open(IN,"<$file[0]") or die "Can't open $file[0] file!\n";
$file[1]=$file[0].".fa";
open(OUT,">$file[1]") or die "Can't create $file[1] file!\n";
while(<IN>)
{
	$line=$_;
	($pos,$len,$status[0],$status[1])=$line=~/pos\:(\d+)\tlength\:(\d+)\t\+\:(\d)\t\-\:(\d)/;
       	$name=$pos."-".$len."-".$status[0]."-".$status[1];
	$primer=substr($seq,$pos,$len);
        print OUT ">$name\n$primer\n";
}
close IN;
close OUT;

if($common_file)
{
	$file[2]=$file[0]."-common_list.txt";
	open(LIST,">$file[2]") or die "Can't create $file[2] file!\n";
	for($i=0;$i<$num_common;$i++)
	{
		print LIST "$list_common[$i]\t$i\n";
	}
	close LIST;

	$file[3]=$file[0]."-common.txt";
	open(COMMON,">$file[3]") or die "Can't create $file[3] file!\n";
}
if($specific_file || $left)
{
	$file[4]=$file[0]."-specific.txt";
	open(SPECIFIC,">$file[4]") or die "Can't create $file[4] file!\n";
}

##run bowtie2
for($i=0;$i<@list_index;$i++)
{
	$file[5]=$file[0]."-".$i.".sam";
	system("$path_bowtie2 -f -N 1 -a --omit-sec-seq --no-unal -p $thread -x $list_index[$i] -U $file[1] -S $file[5]");

	open(BOWTIE,"<$file[5]") or die "Can't open the $file[5] file!\n";
	while(<BOWTIE>)
	{
		$line=$_;
		@array=split("\t",$line);
		if(@array<9)
		{
			next;
		}
		if($array[5]=~/[ID]/)
		{
			next;
		}
		if(length($array[2])<=300)
		{
			$flag=$array[2];
		}
		else
		{
			$flag=substr($array[2],0,300);
		}
		($mismatch)=$line=~/NM:i:(\d+)/;
		($pos,$len,$status[0],$status[1])=$array[0]=~/^(\d+)\-(\d+)\-(\d)\-(\d)/;
	##common
		if($common_file&&(exists $common{$flag})&&($mismatch<=$par_PrimerCommon))
		{
			$begin=0;
			$stop=0;
			if($mismatch>0)
			{
				($MD)=$line=~/MD:Z:([\w\d]+)/;
				($value)=$MD=~/^(\d+)/;
				if($value<5)
				{
					$begin=1;
				}
				($value)=$MD=~/(\d+)$/;
				if($value<5)
				{
					$stop=1;
				}
			}
			$status[2]=0;
			$status[3]=0;

			if(($array[1]-$array[1]%16)/16%2==1) ##minus
			{
				if($status[0]&&($begin==0))
				{
					$status[3]=1;
				}
				if($status[1]&&($stop==0))
				{
					$status[2]=1;
				}
			}
			else ##plus-plus
			{
				if($status[0]&&($stop==0))
				{
					$status[2]=1;
				}
				if($status[1]&&($begin==0))
				{
					$status[3]=1;
				}
			}
			if($status[2]+$status[3]==0)
			{
				next;
			}
		
			print COMMON "$pos\t$len\t$common{$flag}\t$array[3]\t$status[2]\t$status[3]\n";
			next;
		}
		if($mismatch>$par_PrimerSpecific)
		{
			next;
		}
		if(($specific_file&&(exists $specific{$flag}))||$left)
		{
			$begin=0;
                        $stop=0;
                        if($mismatch>0)
                        {
                                ($MD)=$line=~/MD:Z:([\w\d]+)/;
                                ($value)=$MD=~/^(\d+)/;
                                if($value<5)
                                {
                                        $begin=1;
                                }
                                ($value)=$MD=~/(\d+)$/;
                                if($value<5)
                                {
                                        $stop=1;
                                }
                        }
                        $status[2]=0;
                        $status[3]=0;

                        if(($array[1]-$array[1]%16)/16%2==1) ##minus
                        {
                                if($status[0]&&($begin==0))
                                {
                                        $status[3]=1;
                                }
                                if($status[1]&&($stop==0))
                                {
                                        $status[2]=1;
                                }
                        }
                        else ##plus-plus
                        {
                                if($status[0]&&($stop==0))
                                {
                                        $status[2]=1;
                                }
                                if($status[1]&&($begin==0))
                                {
                                        $status[3]=1;
                                }
                        }
                        if($status[2]+$status[3]==0)
                        {
                                next;
                        }
			if(exists $specific{$flag})
			{
				print SPECIFIC "$pos\t$len\t$specific{$flag}\t$array[3]\t$status[2]\t$status[3]\n";
			}
			else
			{
				$specific{$flag}=$num_specific;
				print SPECIFIC "$pos\t$len\t$num_specific\t$array[3]\t$status[2]\t$status[3]\n";
				$num_specific++;
			}
		}
	}
	close BOWTIE;	
	system("rm $file[5]");
}
system("rm $file[1]");
if($common_file)
{
	close COMMON;
}
if($specific_file || $left)
{
	close SPECIFIC;
}
$end=time();
$take=$end-$start;
print "    In this step, it takes $take seconds.\n";

##Probe
if($NoProbe==1)
{
	exit;
}
#####################################run bowtie2 and analysis for probes
print "Now the program is handling the single probe file\n";
$start=time();
$file[0]=$dir."Probe/".$prefix;
open(IN,"<$file[0]") or die "Can't open $file[0] file!\n";
$file[1]=$file[0].".fa";
open(OUT,">$file[1]") or die "Can't create $file[1] file!\n";
while(<IN>)
{
        $line=$_;
        ($pos,$len,$status[0],$status[1])=$line=~/pos\:(\d+)\tlength\:(\d+)\t\+:(\d)\t\-:(\d)/;
	if($status[0]==1)
	{
        	$name=$pos."-".$len."-+";
	}
	else
	{
		$name=$pos."-".$len."--";
	}
        $primer=substr($seq,$pos,$len);
        print OUT ">$name\n$primer\n";
}
close IN;
close OUT;

if($common_file)
{
        $file[3]=$file[0]."-common.txt";
        open(COMMON,">$file[3]") or die "Can't create $file[3] file!\n";
}
if($specific_file || $left)
{
	$file[4]=$file[0]."-specific.txt";
	open(SPECIFIC,">$file[4]") or die "Can't create $file[4] file!\n";
}

##run bowtie2
for($i=0;$i<@list_index;$i++)
{
        $file[5]=$file[0]."-".$i.".sam";
        system("$path_bowtie2 -f -N 1 -a --omit-sec-seq --no-unal -p $thread -x $list_index[$i] -U $file[1] -S $file[5]");

        open(BOWTIE,"<$file[5]") or die "Can't open the $file[5] file!\n";
        while(<BOWTIE>)
        {
                $line=$_;
                @array=split("\t",$line);
                if(@array<9)
                {
                        next;
                }
                if($array[5]=~/[ID]/)
                {
                        next;
                }
                if(length($array[2])<=300)
                {
                        $flag=$array[2];
                }
                else
                {
                        $flag=substr($array[2],0,300);
                }
                ($mismatch)=$line=~/NM:i:(\d+)/;
                ($pos,$len,$status[0])=$array[0]=~/^(\d+)\-(\d+)\-(.)/;

		$status[2]=0;
		$status[3]=0;
	##the T-pos in alignment
		if($status[0] eq "+")
		{
			if(($array[1]-$array[1]%16)/16%2==1) ##minus
			{
				$status[3]=1;
				$begin=15;
				$stop=$len-29;
			}
			else
			{
				$status[2]=1;
				$begin=30;
				$stop=$len-14;
			}
		}
		else
		{
			if(($array[1]-$array[1]%16)/16%2==1) ##minus
			{
				$status[2]=1;        
				$begin=30;
				$stop=$len-14;
			}        
                        else
                        {        
                                $status[3]=1;        
                                $begin=15;
                                $stop=$len-29;
                        }               
                }
        ##common
                if($common_file&&(exists $common{$flag})&&($mismatch<=$par_ProbeCommon))
                {
			($MD)=$line=~/MD:Z:([\w\d]+)/;
			$status[1]=0; ##as a flag
			$MisPos=0;
			while(1)
			{
				if($MD=~/^\d+$/)
				{
					last;
				}
				if($MD=~/^\d+/)
				{
                                	($value)=$MD=~/^(\d+)/;
					$MisPos=$MisPos+$value;
					($MD)=$MD=~/^\d+(.+)$/;
					next;
				}
				$MisPos++;
				if((abs($MisPos-$begin)<5)||(abs($MisPos-$stop)<5))
				{
					$status[1]++;
					last;
				}
				($MD)=$MD=~/^.(.+)$/;
			}
			if($status[1]>0)
			{
				next;
			}
			if($status[2]==1)
			{
				print COMMON "$pos\t$len\t$common{$flag}\t$array[3]\t1\t0\n";
			}
			else
			{
				print COMMON "$pos\t$len\t$common{$flag}\t$array[3]\t0\t1\n";
			}
                        next;
                }

                if($mismatch>$par_ProbeSpecific)
                {
                        next;
                }
                if(($specific_file&&(exists $specific{$flag}))||$left)
                {
			($MD)=$line=~/MD:Z:([\w\d]+)/;
                        $status[1]=0;
                        $MisPos=0;
                        while(1)
                        {
                                if($MD=~/^\d+$/)
                                {
                                        last;
                                }
                                if($MD=~/^\d+/)
                                {
                                        ($value)=$MD=~/^(\d+)/;
                                        $MisPos=$MisPos+$value;
                                        ($MD)=$MD=~/^\d+(.+)$/;
                                        next;
                                }
                                $MisPos++;
                                if((abs($MisPos-$begin)<5)||(abs($MisPos-$stop)<5))
                                {
                                        $status[1]++;
                                        last;
                                }
                                ($MD)=$MD=~/^.(.+)$/;
                        }
                        if($status[1]>0)
                        {
                                next;
                        }
                        if(exists $specific{$flag})
                        {
                                print SPECIFIC "$pos\t$len\t$specific{$flag}\t$array[3]\t";
				if($status[2]==1)
				{
					print SPECIFIC "1\t0\n";
				}
				else
				{
					print SPECIFIC "0\t1\n";
				}
                        }
                        else
                        {
                                $specific{$flag}=$num_specific;
                                print SPECIFIC "$pos\t$len\t$num_specific\t$array[3]\t";
				if($status[2]==1)
                                {
                                        print SPECIFIC "1\t0\n";
                                }
                                else
                                {
                                        print SPECIFIC "0\t1\n";
                                }
                                $num_specific++;
                        }
                }
        }
        close BOWTIE;
        system("rm $file[5]");
}
system("rm $file[1]");
if($common_file)
{
	close COMMON;
}
if($specific_file || $left)
{
	close SPECIFIC;
}
$end=time();
$take=$end-$start;
print "    In this step, it takes $take seconds.\n";
