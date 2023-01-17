#!/usr/bin/perl 
use Math::Trig;
#########################################################################
#	Program: ArgStack_2						#
#									#
#	Jason Thomas Maynes	- Mark Glover				#
#	May/June 2002							#	
#									#	
#	This program will find Arg's that H-bond to a guanine in a 	#
#	pyrimidine-guanine pair and also are within a cutoff distace 	#
#	of the C5 atom of the pyrimidine ring. Will also print out	#
#	statisitics for bond distances that are indicative of a 	#
#	loss of pyrimidine-guanine base stacking.			#
#########################################################################

if ($#ARGV == -1 )
{
	print "\nUsage ./FindArg [file with pdb names - one per line] [H-bond cutoff distance] [Hydrophobic interaction distance]\n\n";
	exit;
}

print "****************************************************************\n\n";
print "				ArgStack_2D				\n\n";
print "		    Jason Thomas Maynes/ Mark Glover - May 2002	      	\n\n";
print "                                                                \n\n";
print "****************************************************************\n\n";

#read in file list
open(PDBFILES,"<$ARGV[0]");
@files=<PDBFILES>;

#open output files
#first, outall contains all data
$output1="outarg";
open (OUTARG,">$output1");

print OUTARG "\nYour H-bonding distance: $ARGV[1]\n";
print OUTARG "Your hydrophobic distance: $ARGV[2]\n\n";
print OUTARG "file	Argchain	residue	CZC5	angle	C2C5	O2C4	O2C5\n";

#process each file
for ($currentpdb=0; $currentpdb <= $#files; $currentpdb++)
{
	open(CURRENTPDBFILE,"<$files[$currentpdb]");
	my @currentpdb=<CURRENTPDBFILE>;
	close(CURRENTPDBFILE);

	#extract residues
	my @residues;
	foreach (@currentpdb)
	{
		my @array=split(" ",$_);

		if ($array[0] eq "ATOM" || $array[0] eq "HETATM")
		{
			push(@residues,$_);
		}
		@array=();
	}

	#look for PyG steps
	my $foundPy=0;	#boolean for if currently found pyrimidine
	my @foundPyG; #array for coord values for current pyrimidine
	my %GHBcoord; #hash for coord values of guanines 
	my %GSTcoord; #hash for coord values of guanines
	my %ARGcoord; #hash for coord values of NH1, NH2, NE of Arg's
	my %ARGCzcoord; #hash for coord values of CZ of Arg's
	my %residuesdone; #hash to count  no. hb's
	my $storePynum; #scalar for current pyrimidine values residue num
	my $storePychain; #scalar for current pyrimidine values residue chain
	my $foundG=0; #boolean for if found a guanine after a pyrimidine
	
	#elementary FSM for finding Py-G sets 
	foreach (@residues)
	{
		my @array=split(" ",$_);

		my $residue = $array[3];
	 	my $atom = $array[2];

		my $chain;
		my $resnum;
		my $coordx;
		my $coordy;
		my $coordz;
		#test for residue #'s over 1000
		{
			my $resnumid1 = $array[5];
			my $resnumid2 = sprintf "%.0f",$resnumid1;
			
			if ($resnumid1 == $resnumid2)
			{
				$chain=$array[4];
				$resnum=$array[5];
				$coordx=$array[6];
				$coordy=$array[7];
				$coordz=$array[8];
			}
			else
			{
				$chain=substr($array[4],0,1);
				$resnum=substr($array[4],1);
				$coordx=$array[5];
				$coordy=$array[6];
				$coordz=$array[7];
			}
		}
		
		#if prev found Py but no G reset
		if (!($residue eq 'G') && !($residue eq "GUA") && $foundPy == 1)
		{
			$foundPy=0;
		}
		#found a Py
		if ($residue eq 'T' || $residue eq 'C' || $residue eq "CYT" || $residue eq "THY")
		{
			$foundPy=1;
			$storePynum=$resnum;
			$storePychain=$chain;
			$foundG=0;
		}
		#if prev found a Py and now a G, store the chain id and residue # of Py and G
		if (($residue eq 'G' || $residue eq "GUA") && $foundPy == 1 && $foundG == 0)
		{
			#ensure that on same chain, also store the coord of centroid of Py ring
			if ($storePynum < $resnum)
			{
				push(@foundPyG,"$chain $storePynum $resnum");
			}
			$foundG=$resnum;
		}
		#get coord of O6 and N7 for guanine and C2/C5 coord for Pyrimidine
		if ($foundG)
		{
			#check for atom, that same chain and that same residue (continuous G's)
			if (($atom eq "N7" || $atom eq "O6") && ($storePynum < $resnum) && $resnum == $foundG)
			{
				$GHBcoord{"$chain $resnum $atom"}="$coordx $coordy $coordz";
			}
		}
		#get coord of C4 and C5 for guanine 
		if ($foundG)
		{
			#check for atom, that same chain and that same residue (continuous G's)
			if (($atom eq "C4" || $atom eq "C5") && ($storePynum < $resnum) && $resnum == $foundG)
			{
				$GSTcoord{"$chain $resnum $atom"}="$coordx $coordy $coordz";
			}
		}
		#get coord of all arginine NH1, NH2 or NZ, initialize residuesdone
		if ($residue eq "ARG")
		{
			if ($atom eq "NH1" || $atom eq "NH2" || $atom eq "NE")
			{
				$ARGcoord{"$chain $resnum $atom"} = "$coordx $coordy $coordz"; 
				$residuesdone{"$chain $resnum"} = 0;
			}
			if ($atom eq "CZ")
			{
				$ARGCzcoord{"$chain $resnum $atom"} = "$coordx $coordy $coordz"; 
			}
		}
		
	}

	#calculate the distances from every N* in Args to O6 and N7 in guanines
	my @closeArgs; #array to hold chain and id of close guanines and Arg's
	my %residuesdone; #hash to check if guanine or Arg has been already included
	foreach $guanine (keys %GHBcoord)
	{
		my @guanarray = split(" ",$GHBcoord{$guanine});
		my $guanx = $guanarray[0];
		my $guany = $guanarray[1];
		my $guanz = $guanarray[2];
		@guanarray = split(" ",$guanine);
		my $guanchain=$guanarray[0];
		my $guannum=$guanarray[1];
		
		foreach $arg (keys %ARGcoord)
		{
			my @argarray = split(" ",$ARGcoord{$arg});
			my $argx=$argarray[0];
			my $argy=$argarray[1];
			my $argz=$argarray[2];
			@argarray=split(" ",$arg);
			my $argchain=$argarray[0];
			my $argnum=$argarray[1];
			
			$distance = sqrt(($guanx-$argx)**2+($guany-$argy)**2+($guanz-$argz)**2);
		
			if ($distance < $ARGV[1] && ($residuesdone{"$argchain $argnum"} == 0))
			{		
				$residuesdone{"$argchain $argnum"}++;
			}
			if ($distance < $ARGV[1] && ($residuesdone{"$argchain $argnum"} == 1))
			{
				push(@closeArgs,"$guanchain $guannum $argchain $argnum");			
				$residuesdone{"$argchain $argnum"}++;
			}
			@argarray=();
		}
		@guanarray=();
	}

	#calculate now if successful Arg's are near Py ring, if near, calculate base stacking distances
	my @successArgs; #Arg's that have been determined to be close to guanine, 5'Py, show stacking dist's
	foreach $checkarg (@closeArgs)
	{
		my @argarray = split (" ",$checkarg);
		my $argchain = $argarray[2];
		my $argnum = $argarray[3];
		my $argx;
		my $argy;
		my $argz;
		my $argxNH1;
		my $argyNH1;
		my $argzNH1;
	
		
		#get coord of CZ for this Arg
		foreach $argcz (keys %ARGCzcoord)
		{
			my @argczarray = split (" ",$argcz);
			if ($argczarray[0] eq $argchain && $argczarray[1] eq $argnum)
			{
				my @argczcoordarray = split (" ",$ARGCzcoord{$argcz});
				$argx=$argczcoordarray[0];
				$argy=$argczcoordarray[1];
				$argz=$argczcoordarray[2];

				my $argNH1coord = $ARGcoord{"$argchain $argnum "."NH1"};
				my @argNH1coordarray=split (" ",$argNH1coord);
				$argxNH1=$argNH1coordarray[0];
				$argyNH1=$argNH1coordarray[1];
				$argzNH1=$argNH1coordarray[2];
			}
		}

		#find C5, C2, O2 of matching pyrimidine, C4/5 of matching G
		my $guanchain = $argarray[0];
		my $guannum = $argarray[1];
		my $pyC5x;
		my $pyC5y;
		my $pyC5z;
		my $pyC2x;
		my $pyC2y;
		my $pyC2z;
		my $pyO2x;
		my $pyO2y;
		my $pyO2z;

		#need to find coord of C5, C2, O2 of Pyrimidine ring
		foreach $LINE (@residues)
		{
			my @array=split (" ",$LINE);
			my $testchain;
			my $testresnum;
			my $coordx;
			my $coordy;
			my $coordz;
	 		my $testatom = $array[2];
			#test for residue #'s over 1000
			{
				my $resnumid1 = $array[5];
				my $resnumid2 = sprintf "%.0f",$resnumid1;
			
				if ($resnumid1 == $resnumid2)
				{
					$testchain=$array[4];
					$testresnum=$array[5];
					$coordx=$array[6];
					$coordy=$array[7];
					$coordz=$array[8];
				}
				else
				{
					$testchain=substr($array[4],0,1);
					$testresnum=substr($array[4],1);
					$coordx=$array[5];
					$coordy=$array[6];
					$coordz=$array[7];
				}
			}
			if ($testchain eq $guanchain && $testresnum == ($guannum-1))
			{
				if ($testatom eq "C5")
				{
					$pyC5x=$coordx;		
					$pyC5y=$coordy;		
					$pyC5z=$coordz;
				}
				if ($testatom eq "C2")
				{
					$pyC2x=$coordx;		
					$pyC2y=$coordy;		
					$pyC2z=$coordz;
				}
				if ($testatom eq "O2")
				{
					$pyO2x=$coordx;		
					$pyO2y=$coordy;		
					$pyO2z=$coordz;
				}
			}
		}
		
		#calculate distance between CZ and C5
		my $distance = sqrt(($pyC5x-$argx)**2+($pyC5y-$argy)**2+($pyC5z-$argz)**2);

		#calculate vector b/w CZ and NH1
		my $argvecx = $argxNH1-$argx;
		my $argvecy = $argyNH1-$argy;
		my $argvecz = $argzNH1-$argz;

		#calculate vector b/w CZ and centroid of ring
		my $pyvecx = $pyC5x-$argx;
		my $pyvecy = $pyC5y-$argy;
		my $pyvecz = $pyC5z-$argz;

		#calculate distance between CZ and NH1
		my $NH1distance = sqrt(($argxNH1-$argx)**2+($argyNH1-$argy)**2+($argzNH1-$argz)**2);

		#calculate the angle between the two vectors
		my $angle = (($argvecx*$pyvecx)+($argvecy*$pyvecy)+($argvecz*$pyvecz))/($NH1distance*$distance);
		$angle = acos($angle);
		$angle = (360*$angle)/(2*3.1415);
		$angle = 90 - $angle;

		#check for cutoff distance and calc stacking distances and save
		if ($distance < $ARGV[2])
		{
		    foreach $guanine (keys %GSTcoord)
		    {
			my @guanarray = split(" ",$GSTcoord{$guanine});
			my $guanx = $guanarray[0];
			my $guany = $guanarray[1];
			my $guanz = $guanarray[2];
			@guanarray = split(" ",$guanine);
			my $testchain=$guanarray[0];
			my $testnum=$guanarray[1];
			my $testatom=$guanarray[2];

		# calc. c2 - c5 distance, o2 - c5 distance
			if ($testchain eq $guanchain && $testnum == ($guannum) && $testatom eq "C5")
			{
			$distc2c5 = sqrt(($pyC2x-$guanx)**2+($pyC2y-$guany)**2+($pyC2z-$guanz)**2);
			$disto2c5 = sqrt(($pyO2x-$guanx)**2+($pyO2y-$guany)**2+($pyO2z-$guanz)**2);
		        }
		# calc. o2 - c4 distance
			if ($testchain eq $guanchain && $testnum == ($guannum) && $testatom eq "C4")
			{
			$disto2c4 = sqrt(($pyO2x-$guanx)**2+($pyO2y-$guany)**2+($pyO2z-$guanz)**2);
		        }
		    }
			push(@successArgs,"$argchain $argnum $distance $angle $distc2c5 $disto2c5 $disto2c4");
		}
	}

	$file = $files[$currentpdb];
	chomp $file;
	#print OUTARG "For File $file: \n";

	#print OUTARG "\tAll Args within H-bonding of guanine in Py-G step:\n";
	foreach (@closeArgs)
	{
		my @array=split(" ",$_);
		#print OUTARG "\t\tArg chain: $array[2] residue: $array[3]\n"
	}
	#print OUTARG "\tAll Args within H-bonding and stacking:\n";
	foreach (@successArgs)
	{
		my @array=split(" ",$_);
		my $distance = sprintf "%.2f", $array[2];
		my $angle = sprintf "%0.2f", $array[3];
		my $d2 = sprintf "%.2f", $array[4];
		my $d3 = sprintf "%.2f", $array[5];
		my $d4 = sprintf "%.2f", $array[6];

		#print OUTARG "\t\tArg chain: $array[0] residue: $array[1] distance CZ to Py C5: $distance angle: $angle\n";
		#print OUTARG "\t\tStacking: C2-C5: $d2 O2-C4: $d4 O2-C5: $d3\n";
		#print OUTARG "file	Argchain	residue	CZC5	angle	C2C5	O2C4	O2C5\n";
		print OUTARG "$file\t $array[0]\t $array[1]\t $distance\t $angle\t $d2\t $d4\t $d3\t\n";
	 }
	#print OUTARG "\n";
	#clear all arrays/hashes for next pdb file
	@currentpdb=();
	@residues=();
	@foundPyG=();
	%GHBcoord=();
	%GSTcoord=();
	%ARGcoord=();
	%ARGCzcoord=();
	@closeArgs=();
	@successArgs=();
	@stack=();
	%residuesdone=();
}

