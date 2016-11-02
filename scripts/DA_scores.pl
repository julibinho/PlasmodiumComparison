#!/usr/bin/perl -w
use strict;
if (exists($ENV{'PERL_LIB'})){ use lib $ENV{'PERL_LIB'};}else{print "ERROR:: Not existis PERL_LIB\n";exit;}


use File::Files;
use Util::Array;
use analisys::analisys;
use Sequence::Domain;


my $PROG_NAME = "Protozoa::scripts::getArchsGroups.pl";

my $verbose    	= "false";
my $archFile  		= "";	
my $groupFile		= "";
my $colGroupName = "";
my $colGroupProtein = "";
my $colArchProtein = "";
my $colArchDomain = "";

my $SEP = "\t";
my %hashGroup	= ();				
my %hashArch  = ();	
my %hashDomains = ();	
my %hashPerc = ();		
my $info = "";	
my $infoGroup = "";	
my $countProteinGroup = 0;
my $countGroupOut = 0;
my $totalGroup = 0;
my $total100Group = 0;
my $outputFile = "";
my $countGroupNotFound = 0;

&getParameters;
&checkParameters;
&main;

############################ sub getParameters ##################################
sub getParameters {
	my $counter;
	my $ending_value = scalar(@ARGV);		#numero de parametros de entrada
	for($counter = 0 ; $counter < $ending_value ; $counter++){
		if($ARGV[$counter] eq "-verbose"){
			$verbose = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-archFile"){
			$archFile = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-groupFile"){
			$groupFile = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-colGroupName"){
			$colGroupName = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-colGroupProtein"){
			$colGroupProtein = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-colArchProtein"){
			$colArchProtein= $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-colArchDomain"){
			$colArchDomain= $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-outputFile"){
			$outputFile= $ARGV[$counter+1];
		}						
	}
}
############################ sub checkParameters ##################################
sub checkParameters {
	print "\n=============================================\n";
	print "Cheking Parameters $PROG_NAME\n";
	
	String->isEmpty($archFile, "ERROR::".$PROG_NAME.":: -archFile is required");
	FilesIO->checkFile($archFile, "ERROR::".$PROG_NAME.":: File NOT found ".$archFile);	
	String->isEmpty($groupFile, "ERROR::".$PROG_NAME.":: Parameter -groupFile is required\n");
	FilesIO->checkFile($groupFile, "ERROR::".$PROG_NAME.":: File NOT found ".$groupFile);	
	String->isEmpty($outputFile, "ERROR::".$PROG_NAME.":: -outputFile is required");

	print "-verbose = $verbose\n";
	print "-archFile = $archFile\n";
	print "-groupFile = $groupFile\n";	
	print "ok\n";	
}

############################ sub main ##################################
sub main{		
			
	my @array_aux = ();				my @archPred;					


    %hashGroup = FilesIO->fileToHash($groupFile, $SEP , $colGroupName);
    %hashArch = FilesIO->fileToHash($archFile, $SEP , $colArchProtein);
    
    while ( my ($key, $value) = each(%hashGroup) ) {
    	if ($key ne "" && $key ne "NONE"){
    		$totalGroup++;
    		#get Domains
    		%hashDomains = ();    		
    		%hashDomains = countDomains($value, $key);
    		%hashPerc = ();	
    		if (keys %hashDomains){
    			%hashPerc = getArchPerc($countProteinGroup);
			}
    		
    		if (keys %hashPerc){
    			$info .= "\n";
    			@array_aux = keys(%hashPerc);
    			@array_aux = sort {$b <=> $a} @array_aux;
    	   		@archPred = $hashPerc{$array_aux[0]};
    	   		@archPred = sort(@archPred);
    	   		if ($array_aux[0] ==  100){
    				$info .= "Consensus Architecture (100%) =>".Array->printArray(",", \@archPred)."\n";
    				$total100Group++;
    			}else{
    				$info .= "Consensus Architecture (100%) => None\n";
    				$info .= "Consensus Architecture (".$array_aux[0].")% =>".Array->printArray(",", \@archPred)."\n";
    			}
    			$info .= "//\n";
    		}
    	}    	
    }
    
    FilesIO->saveFile($outputFile, $info, $verbose);
    FilesIO->saveFile($outputFile.".groups", $infoGroup, $verbose);
   
    print "Total Groups = $totalGroup\n";
    print "Total Groups Not found = $countGroupOut\n";
    print "Total Groups Never found = $countGroupNotFound\n";
    print "Total Groups agree 100% = $total100Group\n";

}


############################ sub main ##################################
sub calculeDAscore{	
	my $proteinsGroup = $_[0];
	my $groupName = $_[1];
	
	print $groupName."\n";
	
	my $protein;				my @arrayDomains = ();
	my $countNotFound = 0;		$countProteinGroup = 0;
	my $infoProt = "";
	
	my @array_aux = String->StringToArray2d($proteinsGroup, "\n", $SEP);
	$infoProt .= "Group=$groupName\n";
	for (my $j =0; $j < scalar(@array_aux); $j++) { 
    	$protein = $array_aux[$j][$colGroupProtein];  
    	$countProteinGroup++;  		
    	#print $protein."\n";
    	if (exists($hashArch{$protein})){
    		@arrayDomains = String->StringToArray2d($hashArch{$protein}, "\n", $SEP);
    		@arrayDomains = Array->copyColMatrixToArray2($colArchDomain, \@arrayDomains);
    		@arrayDomains = sort(@arrayDomains);
    		$infoProt .= "$protein -". Array->printArray(",", \@arrayDomains)."\n";
    		@arrayDomains = Array->removeRepeated(\@arrayDomains);
    		for (my $k =0; $k < scalar(@arrayDomains); $k++) { 
    			if (exists($hashDomains{$arrayDomains[$k]})){
    				$hashDomains{$arrayDomains[$k]} = $hashDomains{$arrayDomains[$k]} + 1;
    			}else{
    				$hashDomains{$arrayDomains[$k]} = 1;
    			}
    		}
    	}else{
    		$countNotFound++;
    	}
    }
    if ($countNotFound == $countProteinGroup){
    	$countGroupNotFound++;
    }
    $infoGroup .= $groupName."\t".$countNotFound."\t". $countProteinGroup."\n";
    if ($countNotFound == 0){
    	$info .= $infoProt;
    }else{
    	%hashDomains = ();
    	$countGroupOut++;
    }
    return %hashDomains;
    
}


############################ sub main ##################################
sub calculeDAscore{	
	my $TotalProtein= $_[0];
	my $perc;
	
	$info .= "Group Archs =>";
    my @archPred = (); 
    while ( my ($keyDomain, $valueDomain) = each(%hashDomains) ) {
    	$perc = $valueDomain/$TotalProtein;
    	$perc = sprintf("%.2f", $perc)*100;
    	$info .=  $keyDomain."(".$perc. "), ";
    	if (exists($hashPerc{$perc})){
    		$hashPerc{$perc} = $hashPerc{$perc}.", ".$keyDomain;
    	}else{
    		$hashPerc{$perc} = $keyDomain;
    	}    		
    }
	return %hashPerc;
}