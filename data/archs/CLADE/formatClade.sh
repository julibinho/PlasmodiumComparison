for ORG in yoelii17X berghei yoeliiYM reichenowiCDC knowlesiH falciparumIT cynomolgi chabaudi yoelii17XNL vivax;
do
	echo $ORG;
    perl ~/projects/perlLib/run/getColunsFile.pl -file p.$ORG.txt -cols 5:1:2:0:4:7 -outputFile p.$ORG.clade
    
done

