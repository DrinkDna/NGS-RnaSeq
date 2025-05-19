#!/usr/bin/bash

#________________________________________________________________________________________________________________________________________________

echo "quality checking"

#before trim 

mkdir -p QC
#fastqc ./*.fastq.gz -o ./QC

#after trim 
mkdir -p trimmed_qc

fastqc ./trimmed/*_paired.fq.gz -o ./trimmed_qc

#_________________________________________________________________________________________________________________________________________________

echo "quality trimming..."

#!/usr/bin/bash

mkdir -p trimmed

for file1 in *_1.fastq.gz; do
    # Extract base filename (without _1.fastq.gz)
    base=$(basename "$file1" _1.fastq.gz)
    file2="${base}_2.fastq.gz"

    # Define output file names
    out1p="trimmed/${base}_1_paired.fq.gz"
    out1u="trimmed/${base}_1_unpaired.fq.gz"
    out2p="trimmed/${base}_2_paired.fq.gz"
    out2u="trimmed/${base}_2_unpaired.fq.gz"

    # Run Trimmomatic
    java -jar trimmomatic-0.39.jar PE -phred33 -threads 4 \
        "$file1" "$file2" \
        "$out1p" "$out1u" \
        "$out2p" "$out2u" \
        ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 MINLEN:36
done


#_____________________________________________________________________________________________________________________________________________________

#!/usr/bin/bash

echo "build genome index"

index="./index"
mkdir -p "$index"

# Build HISAT2 index
hisat2-build /home/jsspsls/suraj/sbmtech/PRJ2025JOB/GCF_000001405.40_GRCh38.p14_genomic.fna /home/jsspsls/suraj/sbmtech/PRJ2025JOB/index/grch38index

#________________________________________________________________________________________________________________________________________________________

#!/usr/bin/bash

echo "Alignment with reference genome using custom scoring penalties"

trimmed="/home/jsspsls/suraj/sbmtech/PRJ2025JOB/trimmed"
aligned="/home/jsspsls/suraj/sbmtech/PRJ2025JOB/aligned"
hisat2_index="/home/jsspsls/suraj/sbmtech/PRJ2025JOB/index/grch38index"

mkdir -p "$aligned"

for file in "$trimmed"/*_1_paired.fq.gz; do
    base_name=$(basename "$file" "_1_paired.fq.gz")

    file_1="$trimmed/${base_name}_1_paired.fq.gz"
    file_2="$trimmed/${base_name}_2_paired.fq.gz"
    log="$aligned/${base_name}.log"

    echo "[$(date)] Aligning $base_name with custom scoring..." | tee "$log"

    hisat2 -x "$hisat2_index" \
           -1 "$file_1" -2 "$file_2" \
           --mp 6,2 \
           --np 1 \
           --rdg 5,3 \
           --score-min L,0,-0.2 \
           2>>"$log" \
    | samtools view -Sb - 2>>"$log" \
    | samtools sort -o "$aligned/${base_name}_sorted.bam" 2>>"$log"
    
    samtools index "$aligned/${base_name}_sorted.bam" 2>>"$log"

    echo "[$(date)] $base_name alignment done." | tee -a "$log"
done

echo "All alignments with custom scoring completed."

#____________________________________________________________________________________________________________________________________________________________


#!/usr/bin/bash

#infer experiment

gtf2bed < /media/shubham/NV21/PRJ2025JOB/GCF_000001405.40_GRCh38.p14_genomic.gtf > /media/shubham/NV21/PRJ2025JOB/aligned/annotation.bed
infer_experiment.py -i /media/shubham/NV21/PRJ2025JOB/aligned/SRR33334923_sorted.bam -r /media/shubham/NV21/PRJ2025JOB/aligned/annotation.bed

#______________________________________________________________________________________________________________________________________________________________

#!/usr/bin/bash

echo "extracting counts"

mkdir -p counts

featureCounts -T 8 \
    -a ./GCF_000001405.40_GRCh38.p14_genomic.gtf \
    -o ./countsunstrand/fcounts.txt \
    -s 0 \
    -p \
    --countReadPairs \
    /media/shubham/NV21/PRJ2025JOB/aligned/*_sorted.bam
    
#________________________________________________________________________________________________________________________________________________________________

#!/usr/bin/bash

echo "merge counts"

# Extract header (sample names) from fcounts.txt
head -n 1 fcounts.txt > output.csv

# Extract only gene_id and counts (from column 7 onward)
awk 'NR>1 { 
    printf "%s", $1;
    for (i=7; i<=NF; i++) {
        printf ",%s", $i
    }
    print ""
}' fcounts.txt >> output.csv


#__________________________________________________________________________________________________________________________________________________________________

#DESeq2 Analysis in R with count matrix and sample information 
























