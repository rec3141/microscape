#!/usr/bin/env bash

# input DIRECTORY as $1
# must have cutadapt installed in $PATH

# takes input after initial demultiplexing of outer barcode on MiSeq

for F1 in `find $1 -name "*fastq.gz" | cut -f1-2 -d'_' | uniq`; do 
echo $F1

DIR=`dirname $F1`
FILE=`basename $F1`

for FQ in `find $DIR -name "$FILE*R1*.fastq.gz" | grep -v trimmed`; do
	RQ=${FQ/_R1/_R2}

	TFO=${FQ/fastq\.gz/ctmp\.fastq\.gz}
	TRO=${RQ/fastq\.gz/ctmp\.fastq\.gz}
	FQO=${FQ/fastq\.gz/ctrimmed\.fastq\.gz}
	RQO=${RQ/fastq\.gz/ctrimmed\.fastq\.gz}
	ST=`basename $FQO .fastq.gz`

	PRIMERS=primers-all.fa;

	if [[ "$ST" =~ ^16.* ]]; 
	    then PRIMERS=primers-bac.fa;
	    elif [[ "$ST" =~ ^18.* ]]; 
	    then PRIMERS=primers-euk.fa; 
	    elif [[ "$ST" =~ ^ITS.* ]];
	    then PRIMERS=primers-its.fa;
	    else PRIMERS=primers-all.fa; 
	fi

    cutadapt -g file:$PRIMERS -G file:$PRIMERS --discard-untrimmed -j 8 -e 0.12 -o $FQO -p $RQO $FQ $RQ

	done;

	rm $TFO $TRO
done;
