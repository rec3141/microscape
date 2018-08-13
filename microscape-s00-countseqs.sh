# for making plate heatmaps for QC

rm count-demux.txt

for F1 in `ls *.fastq.gz | cut -f1-2 -d'_' | sort -u`; do
 echo "counting $F1";
 A1=`gunzip -c $F1*R1.fastq.gz | wc -l`
 A2=`gunzip -c $F1*R2.fastq.gz | wc -l`
printf '%s\t%s\t%s\n' "$F1" "$A1" "$A2"  >> count-demux.txt 
done;

Rscript microscape-s00-plot-seq-map.R
