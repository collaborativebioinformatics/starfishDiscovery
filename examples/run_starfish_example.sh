#!/usr/bin/env bash

set -euo pipefail

# Run the starfish example

STARFISHDIR=/opt/starfish/

cp -r examples/ test/
cd test/

realpath assembly/* | perl -pe 's/^(.+?([^\/]+?).fasta)$/\2\t\1/' > ome2assembly.txt
realpath gff3/* | perl -pe 's/^(.+?([^\/]+?).final.gff3)$/\2\t\1/' > ome2gff.txt

cat gff3/*.gff3 > macpha6.gff3

if [ ! -d "blastdb" ]; then
  mkdir blastdb
fi
cut -f2 ome2assembly.txt | xargs cat > blastdb/macpha6.assemblies.fna
makeblastdb -in blastdb/macpha6.assemblies.fna -out blastdb/macpha6.assemblies -parse_seqids -dbtype nucl

$STARFISHDIR/aux/seq-gc.sh -Nbw 1000 blastdb/macpha6.assemblies.fna > macpha6.assemblies.gcContent_w1000.bed
rm blastdb/macpha6.assemblies.fna

cut -f1,12  ann/*emapper.annotations | grep -v  '#' | grep -v -P '\t-' | perl -pe 's/\t/\tEMAP\t/' | grep -vP '\tNA' > ann/macph6.gene2emap.txt

cut -f1,10 ann/*emapper.annotations | grep -v '#' | perl -pe 's/^([^\s]+?)\t([^\|]+).+$/\1\t\2/' > ann/macph6.gene2og.txt

$STARFISHDIR/aux/geneOG2mclFormat.pl -i ann/macph6.gene2og.txt -o ann/

if [ ! -d "geneFinder" ]; then
  mkdir geneFinder
fi
echo "According to the starfish example, you should see:"
echo "found 1 new tyr genes and 9 tyr genes that overlap with 11 existing genes"

echo "Running starfish annotate. This may take a while (10 - 20 minutes)"
starfish annotate -T 2 -x macpha6_tyr -a ome2assembly.txt -g ome2gff.txt -p $STARFISHDIR/db/YRsuperfams.p1-512.hmm -P $STARFISHDIR/db/YRsuperfamRefs.faa -i tyr -o geneFinder/

starfish consolidate -o ./ -g macpha6.gff3 -G geneFinder/macpha6_tyr.filt_intersect.gff

realpath macpha6_tyr.filt_intersect.consolidated.gff | perl -pe 's/^/macpha6\t/' > ome2consolidatedGFF.txt

echo "According to the starfish example, you should see:"
echo "found 10 neighborhoods containing input query genes"

echo "Running starfish sketch"

starfish sketch -m 10000 -q geneFinder/macpha6_tyr.filt_intersect.ids -g ome2consolidatedGFF.txt -i s -x macpha6 -o geneFinder/


grep -P '\ttyr\t' geneFinder/macpha6.bed > geneFinder/macpha6.tyr.bed 

# mkdir elementFinder
if [ ! -d "elementFinder" ]; then
  mkdir elementFinder
fi

echo "According to the starfish example, you should see:"
echo "found element boundaries and insertion sites for 7 tyr captains out of 10 input captains"

echo "Running starfish insert"
starfish insert -T 2 -a ome2assembly.txt -d blastdb/macpha6.assemblies -b geneFinder/macpha6.tyr.bed -i tyr -x macpha6 -o elementFinder/


echo "According to the starfish example, you should see:"
echo "found 2 captains with DR boundaries and 5 captains with DR-TIR boundaries out of 7 input captains with candidate DRs"

echo "Running starfish flank"
starfish flank -a ome2assembly.txt -b elementFinder/macpha6.insert.bed -x macpha6 -o elementFinder/

starfish summarize -a ome2assembly.txt -b elementFinder/macpha6.flank.bed -x macpha6 -o elementFinder/ -S elementFinder/macpha6.insert.stats -f elementFinder/macpha6.flank.singleDR.stats -g ome2consolidatedGFF.txt -A ann/macph6.gene2emap.txt -t geneFinder/macpha6_tyr.filt_intersect.ids 

awk '{
    if ($4 ~ /DR/) print $0 "\t255,255,0";        # Yellow for 'DR' in column 4
    else if ($4 ~ /TIR/) print $0 "\t255,165,0";  # Orange for 'TIR' in column 4
    else if ($5 ~ /tyr|cap/) print $0 "\t255,0,0"; # Red for 'tyr' or 'cap' in column 5
    else if ($5 != ".") print $0 "\t128,0,128";   # Purple if column 5 is not '.'
    else print $0 "\t169,169,169";                # Dark gray otherwise
}' elementFinder/macpha6.elements.bed > elementFinder/macpha6.elements.color.bed

hmmsearch --noali --notextw -E 0.001 --max --cpu 12 --tblout elementFinder/macpha6_tyr_vs_YRsuperfams.out $STARFISHDIR/db/YRsuperfams.p1-512.hmm geneFinder/macpha6_tyr.filt_intersect.fas
perl -p -e 's/ +/\t/g' elementFinder/macpha6_tyr_vs_YRsuperfams.out | cut -f1,3,5 | grep -v '#' | sort -k3,3g | awk '!x[$1]++' > elementFinder/macpha6_tyr_vs_YRsuperfams_besthits.txt

grep -P '\tcap\t' elementFinder/macpha6.elements.bed | cut -f4,7 > elementFinder/macpha6.cap2ship.txt
$STARFISHDIR/aux/searchReplace.pl --strict -i elementFinder/macpha6_tyr_vs_YRsuperfams_besthits.txt -r elementFinder/macpha6.cap2ship.txt > elementFinder/macpha6_elements_vs_YRsuperfams_besthits.txt

mmseqs easy-cluster geneFinder/macpha6_tyr.filt_intersect.fas elementFinder/macpha6_tyr elementFinder/ --threads 2 --min-seq-id 0.5 -c 0.25 --alignment-mode 3 --cov-mode 0 --cluster-reassign
$STARFISHDIR/aux/mmseqs2mclFormat.pl -i elementFinder/macpha6_tyr_cluster.tsv -g navis -o elementFinder/

starfish sim -m element -t nucl -b elementFinder/macpha6.elements.bed -x macpha6 -o elementFinder/ -a ome2assembly.txt
starfish group -m mcl -s elementFinder/macpha6.element.nucl.sim -i hap -o elementFinder/ -t 0.05

$STARFISHDIR/aux/searchReplace.pl -i elementFinder/macpha6_tyr_cluster.mcl -r elementFinder/macpha6.cap2ship.txt > elementFinder/macpha6.element_cluster.mcl

$STARFISHDIR/aux/mergeGroupfiles.pl -t elementFinder/macpha6.element_cluster.mcl -q elementFinder/macpha6.element.nucl.I1.5.mcl > elementFinder/macpha6.element.navis-hap.mcl

awk '{ for (i = 2; i <= NF; i++) print $i"\t"$1 }' elementFinder/macpha6.element.navis-hap.mcl > elementFinder/macpha6.element.navis-hap.txt

join -t$'\t' -1 1 -2 2 <(sort -t$'\t' -k1,1 elementFinder/macpha6.element.navis-hap.txt | grep -P '_e|_s') <(sort -t$'\t' -k2,2 elementFinder/macpha6.elements.feat) | awk -F'\t' '{print}' > elementFinder/macpha6.elements.temp.feat
echo -e "#elementID\tfamilyID\tnavisHapID\tcontigID\tcaptainID\telementBegin\telementEnd\telementLength\tstrand\tboundaryType\temptySiteID\temptyContig\temptyBegin\temptyEnd\temptySeq\tupDR\tdownDR\tDRedit\tupTIR\tdownTIR\tTIRedit\tnestedInside\tcontainNested" > elementFinder/macpha6.elements.ann.feat
join -t$'\t' -1 1 -2 1 <(sort -t$'\t' -k1,1 elementFinder/macpha6_elements_vs_YRsuperfams_besthits.txt | grep -P '_e|_s' | cut -f1,2) <(sort -t$'\t' -k1,1 elementFinder/macpha6.elements.temp.feat) | awk -F'\t' '{print}' >> elementFinder/macpha6.elements.ann.feat

# mkdir pairViz
if [ ! -d "pairViz" ]; then
  mkdir pairViz
fi
starfish pair-viz -m all -t empty -T 2 -A nucmer -a ome2assembly.txt -b elementFinder/macpha6.elements.bed -f elementFinder/macpha6.flank.singleDR.stats -S elementFinder/macpha6.elements.named.stats -o pairViz/

# mkdir regionFinder
if [ ! -d "regionFinder" ]; then
  mkdir regionFinder
fi

grep -f <(comm -23 <(cut -f1 geneFinder/macpha6_tyr.filt_intersect.ids | sort) <(grep -P '\tcap\t|\ttyr\t' elementFinder/macpha6.elements.bed | cut -f4| sort)) geneFinder/macpha6.tyr.bed > regionFinder/unaffiliated_tyrs.bed

$STARFISHDIR/aux/filterOG.pl -O ann/macph6.gene2og.mcl -a 1 -c 5 -o ann/

echo "According to the starfish example, you should see:"
echo "found 5 regions with at least 1 cross-referenced element-insertion site pair"

echo "Running starfish dereplicate"
starfish dereplicate -e elementFinder/macpha6.element.navis-hap.mcl -t regionFinder/unaffiliated_tyrs.bed -F elementFinder/macpha6.elements.feat -S elementFinder/macpha6.elements.named.stats -O ann/macph6.gene2og.a1.c5.txt -g ome2gff.txt -x macpha6 -o regionFinder/ --flanking 3 --mismatching 1

grep -v '#' regionFinder/macpha6.fog3.d600000.m1.dereplicated.txt | cut -f2 | sort | uniq -c | perl -pe 's/ +//' | sort -k1,1nr

# mkdir locusViz
if [ ! -d "locusViz" ]; then
  mkdir locusViz
fi
starfish locus-viz -T 2 -m region -a ome2assembly.txt -b elementFinder/macpha6.elements.bed -x macpha6 -o locusViz/ -A nucmer -r regionFinder/macpha6.fog3.d600000.m1.regions.txt -d regionFinder/macpha6.fog3.d600000.m1.dereplicated.txt -j regionFinder/macpha6.fog3.d600000.m1.haplotype_jaccard.sim  -g ome2consolidatedGFF.txt --tags geneFinder/macpha6_tyr.filt_intersect.ids --gc macpha6.assemblies.gcContent_w1000.bed

cp -r * /tmp/starfish_example_output/