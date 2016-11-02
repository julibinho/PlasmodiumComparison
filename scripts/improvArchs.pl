#!/usr/bin/perl -w
use strict;
if (exists($ENV{'PERL_LIB'})){ use lib $ENV{'PERL_LIB'};}else{print "ERROR:: Not existis PERL_LIB\n";exit;}


use File::Files;
use Util::Array;
use Align::Primary::blast;
use analisys::analisys;
use Sequence::Domain;


my $PROG_NAME = "Plasmobase::scripts::improvArchs.pl";

my $verbose    		= "false";
my $archFile  		= "";	
my $missingFile		= "";
my $potDomainsFile	= "";
my $colProteinPotDomain	= 3;
my $colProteinArch = 3;
my $colProtMissing = 0;
my $colArchProtein 	= "";
my $colArchDomain 	= "";

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
		}elsif($ARGV[$counter] eq "-missingFile"){
			$missingFile = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-potDomainsFile"){
			$potDomainsFile = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-colProteinPotDomain"){
			$colProteinPotDomain = $ARGV[$counter+1];
		}elsif($ARGV[$counter] eq "-colProteinArch"){
			$colProteinArch= $ARGV[$counter+1];
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
	
	String->isEmpty($missingFile, "ERROR::".$PROG_NAME.":: Parameter -missingFile is required\n");
	FilesIO->checkFile($missingFile, "ERROR::".$PROG_NAME.":: File NOT found ".$missingFile);
		
	String->isEmpty($potDomainsFile, "ERROR::".$PROG_NAME.":: Parameter -potDomainsFile is required\n");
	FilesIO->checkFile($potDomainsFile, "ERROR::".$PROG_NAME.":: File NOT found ".$potDomainsFile);
	
	String->isEmpty($outputFile, "ERROR::".$PROG_NAME.":: -outputFile is required");

	print "-verbose = $verbose\n";
	print "-archFile = $archFile\n";
	print "-missingFile = $missingFile\n";
	print "-potDomainsFile = $potDomainsFile\n";	
	print "-outputFile = $outputFile\n";
	
	print "ok\n";	
}

############################ sub main ##################################
sub main{		
			
	my @array_aux = ();				my @array_missing  = ();					
	my %hashPotDomains;				my %hashArch;
	my $line;						my $codProt;
	my @array_arch = ();	
	my $overlapAA_alowed = 30;		my $perc_overlap_alowed = 0.5;
	my $domains_overlap= "";		my $startHit_k;
	my $endHit_k;					my $domain_k;
	my $startHit_l;					my $endHit_l;					
	my $domain_l;					my $counNoOverlp;
	my $domain;						my $infoDom = "";
	my $info = "";					my $countNewDomains = 0;
	my $infoNews = "";				my $indexArch;
	my $infoPotDomain;
	
	my $colEvA = 0;	my $colStarA = 1;my $colEndA = 2;my $colProtA = 3; my $colDomA = 4; my $colModelA = 5;
    my $colEvPD = 5; my $colStarPD = 1;my $colEndPD = 2;my $colProtPD = 3; my $colDomPD = 4; my $colModelPD = 0;
    
    @array_missing = FilesIO->fileToArray2d($missingFile, "\t");
    %hashPotDomains = FilesIO->fileToHash($potDomainsFile, $SEP , $colProteinPotDomain);
    %hashArch = FilesIO->fileToHash($archFile, $SEP , $colProteinArch);
    
    
    for (my $i = 0; $i < scalar(@array_missing); $i++) {
    	$codProt = $array_missing[$i][$colProtMissing];
    	@array_arch = (); $infoDom = ""; $indexArch = 0;
    	if (exists($hashArch{$codProt})){
    		@array_arch = String->StringToArray2d($hashArch{$codProt}, "\n", "\t");
    		$indexArch = scalar(@array_arch);
    	}
    	
    	for (my $j = 2; $j <= $#{$array_missing[$i]} ; $j++) {
    		#$infoPotDomain = FilesIO->getLinesInFile($potDomainsFile, $codProt, $verbose);
    		
    		#if ($infoPotDomain ne ""){
    		if (exists($hashPotDomains{$codProt})){
    			$infoPotDomain = $hashPotDomains{$codProt};
    			$domain = $array_missing[$i][$j];
    			@array_aux = String->StringToArray2d($infoPotDomain, "\n", "\t");
    			
    			@array_aux = Array->sortArray2d($colEvPD, \@array_aux);
    			 
    						
    			for (my $k = 0; $k < scalar(@array_aux); $k++) {
    				$startHit_k = $array_aux[$k][$colStarPD];  $endHit_k = $array_aux[$k][$colEndPD]; $domain_k = $array_aux[$k][$colDomPD];
    				if ($domain_k eq $domain){
    					$counNoOverlp = 0;
    					for (my $l = 0; $l < scalar(@array_arch); $l++) {
    						$startHit_l = $array_arch[$l][$colStarA];  $endHit_l = $array_arch[$l][$colEndA]; $domain_l = $array_arch[$l][$colDomA];
    						if (Blast->isOverlapHit($startHit_k, $endHit_k, $startHit_l, $endHit_l, $overlapAA_alowed, $perc_overlap_alowed, $domains_overlap, $domain_k, $domain_l) eq "true")	{
    							$l = scalar(@array_arch);
    						}else{
    							$counNoOverlp++;
    						}
    					}
    					if ($counNoOverlp == scalar(@array_arch)){
    						#print "Improve ". $array_aux[$k][$colEvPD]."\t".$array_aux[$k][$colStarPD]."\t".$array_aux[$k][$colEndPD]."\t".$array_aux[$k][$colProtPD]."\t".$array_aux[$k][$colDomPD]."\t".$array_aux[$k][$colModelPD]."\n";
    						$infoDom .= $array_aux[$k][$colEvPD]."\t".$array_aux[$k][$colStarPD]."\t".$array_aux[$k][$colEndPD]."\t".$array_aux[$k][$colProtPD]."\t".$array_aux[$k][$colDomPD]."\t".$array_aux[$k][$colModelPD]."\n";
    						$array_arch[$indexArch][$colEvA] = $array_aux[$k][$colEvPD];
    						$array_arch[$indexArch][$colStarA] = $array_aux[$k][$colStarPD];
    						$array_arch[$indexArch][$colEndA] = $array_aux[$k][$colEndPD];
    						$array_arch[$indexArch][$colProtA] = $array_aux[$k][$colProtPD];
    						$array_arch[$indexArch][$colDomA] = $array_aux[$k][$colDomPD];
    						$array_arch[$indexArch][$colModelA] = $array_aux[$k][$colModelPD];
    						$indexArch++;
    						#print Array->printMatrix2("\t", "\n", \@array_aux);
    						
    						$k = scalar(@array_aux);
    						$countNewDomains++;
    					}
    				}
    			}
    		}
    	}
    	$info .= Array->printMatrix2("\t", "\n", \@array_arch);
    	if ($infoDom ne ""){
       		$infoNews .=$infoDom;
    	}
    }
        
    FilesIO->saveFile($outputFile, $info, $verbose);
    FilesIO->saveFile($outputFile.".news", $infoNews, $verbose);
    print "Add New domains $countNewDomains\n";
    #FilesIO->saveFile($outputFile.".groups", $infoGroup, $verbose);
   
   

}
