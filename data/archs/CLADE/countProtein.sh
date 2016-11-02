for ORG in falciparum yoelii17X berghei yoeliiYM reichenowiCDC knowlesiH falciparumIT cynomolgi chabaudi yoelii17XNL vivax;
do
   echo $ORG;
   cut -f 1 p.$ORG.txt | sort | uniq | wc -l    
done


