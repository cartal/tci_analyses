#!/bin/bash -x

###DATA
REF=/media/prussia/data/refseq/tci_sylvio_x10/sequence/bwaindex/TcI_Sylvio_X10.v2.fasta
GATK_REF=/media/prussia/data/refseq/tci_sylvio_x10/sequence/bwaindex/TcI_Sylvio_X10.v2.fasta
HASH=/media/prussia/data/refseq/tci_sylvio_x10/sequence/stampyindex/TcI_Sylvio_X10
DATA_DIR=

###TOOLS

CPU=
SS=<TOOLS_PATH>/speedseq
BEAGLE=<TOOLS_PATH>/beagle.r1399.jar
PICARD=<TOOLS_PATH>/picard-tools-1.137/picard.jar
STAMPY=<TOOLS_PATH>/stampy-1.0.23/stampy.py
DELLY=<TOOLS_PATH>/delly2/src/delly
TD=<TOOLS_PATH>/TIDDIT/bin/TIDDIT
GATK=<TOOLS_PATH>/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar
FK=<TOOLS_PATH>/fermikit/fermi.kit

for i in $(cat sample_list.txt)

do
    
    mkdir $i
    cd $i

    ###Read filtering 

    cat $DATA_DIR/$i/02-FASTQ/*/*_R1_*.fastq.gz > pair_1.fastq.gz
    cat $DATA_DIR/$i/02-FASTQ/*/*_R2_*.fastq.gz > pair_2.fastq.gz

    nesoni clip --quality 20 --length 45 --gzip yes --homopolymers yes --out-separate yes $i.g.ipe.Q20L45H pairs: pair_1.fastq.gz pair_2.fastq.gz 

    PE1=$i.g.ipe.Q20L45H_R1.fq.gz
    PE2=$i.g.ipe.Q20L45H_R2.fq.gz
    rm pair_1.fastq.gz pair_2.fastq.gz *_single.fq.gz

    ###Primary mapping with BWA

    bwa mem -t $CPU -M $REF $PE1 $PE2 | samtools view -@ $CPU -Sb -> $i.bwa.raw.bam

    ###Secondary mapping with Stampy

    python $STAMPY -g $HASH -h $HASH -t $CPU --bamkeepgoodreads --substitutionrate=0.015 -M $i.bwa.raw.bam | samtools view -@ $CPU -Sb -> $i.stampy.raw.bam

    ###Sort BAM file

    java -jar -Xmx10g $PICARD SortSam  I=$i.stampy.raw.bam SORT_ORDER=coordinate O=$i.sorted.bam

    ###Remove PCR duplicates

    java -jar -Xmx10g $PICARD MarkDuplicates  I=$i.sorted.bam METRICS_FILE=$i-metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT O=$i.s.d.bam

    ###Add RG

    java -jar -Xmx10g $PICARD AddOrReplaceReadGroups I=$i.s.d.bam SORT_ORDER=coordinate RGID=$i RGLB=$i RGPU=SLL-STO RGSM=$i RGPL=Illumina CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT O=$i.RG.bam

    ###Call SNV using GATK-HC

    java -jar -Xmx40g $GATK -T HaplotypeCaller -R $GATK_REF -I $i.RG.bam -nct 24 -o $i.GATK-HC.RAW.vcf

    java -jar -Xmx2g $BEAGLE gtgl=$i.GATK-HC.RAW.vcf nthreads=$CPU impute=true usephase=true out=$i.GATK-HC.RAW.phased

    ###Call SNV using SpeedSeq

    $SS var -t $CPU -o $i.FreeBayes.RAW -q 30 $REF $i.RG.bam

    dtrx $i.FreeBayes.RAW.vcf.gz

    java -jar -Xmx2g $BEAGLE gtgl=$i.FreeBayes.RAW.vcf nthreads=$CPU impute=true usephase=true out=$i.FreeBayes.RAW.phased

    ###Call variants with FermiKit

    $FK/fermi2.pl unitig -s 55m -t $CPU -l 100 -p $i  "$FK/seqtk mergepe $PE1 $PE2 | $FK/trimadap-mt -p4" > $i.mak

    make -f $i.mak -j $CPU

    $FK/run-calling -t $CPU $REF $i.mag.gz | sh

    ###Call SV with LUMPY

    $SS sv -t $CPU -R $REF -o $i.LUMPY.RAW -B $i.RG.bam

    ###Call structural variants with DELLY

    export OMP_NUM_THREADS=$CPU

    samtools index $i.RG.bam

    $DELLY --genome $REF --map-qual 30 --type DEL --outfile $i.DEL.Q25.vcf.gz $i.RG.bam

    $DELLY --genome $REF --map-qual 30 --type INS --outfile $i.INS.Q25.vcf.gz $i.RG.bam

    $DELLY --genome $REF --map-qual 30 --type DUP --outfile $i.DUP.Q25.vcf.gz $i.RG.bam

    $DELLY --genome $REF --map-qual 30 --type INV --outfile $i.INV.Q25.vcf.gz $i.RG.bam

    $DELLY --genome $REF --map-qual 30 --type TRA --outfile $i.TRA.Q25.vcf.gz $i.RG.bam


    ###Call structural variants with FT

    $TD --sv --orientation innie --insert 650 -p 10 -q 30 -b $i.RG.bam -o $i-TIDDIT

    $TD --cov --bin_size 1000 -b $i.RG.bam -o $i-TIDDIT_COV

    cd ..

done
