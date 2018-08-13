# DEMULTIPLEX USING MR_DEMUXY

find . -name "*.gz" -exec mv {} ./ \;
rmdir *
ls *.gz | parallel gunzip {}

for F1 in `ls *.fastq | cut -f1 -d'_' | sort -u`; do 
 echo "demultiplexing $F1"; 
 /work/cryomics/apps/Mr_Demuxy-1.2.0/bin/pe_demuxer.py -r1 $F1*_R1_001.fastq -r2 $F1*_R2_001.fastq -r1_bc forward_bcs.txt -r2_bc reverse_bcs.txt

 cd pe_demuxer_output
 for D2 in R1 R2; do
   cd $D2
   for FILE in *.fastq; do
    mv $FILE $F1"_"$FILE
   done;
   cd ..
 done;

 cd ..

 ls pe_demuxer_output/R1/*.fastq | parallel gzip {}
 ls pe_demuxer_output/R2/*.fastq | parallel gzip {}
 
 mv pe_demuxer_output/*/*.gz ./

 rm -rf pe_demuxer_output

done;