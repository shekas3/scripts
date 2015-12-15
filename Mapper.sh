reference="/home/shahrokh/bio720/FinalProject/reference/Esals173"
annotation="/home/shahrokh/bio720/FinalProject/reference/Esals173.gff3"
samples[1]="SDR_GCCAAT_L002_f"
samples[2]="SD_TGACCA_L002_f"
samples[3]="SRW-1_ACAGTG_L102_f"
samples[4]="SRW-2_CAGATC_L102_f"
samples[5]="SWW-3_CTTGTA_L102_f"
samples[6]="SWW_CTTGTA_L002_f"

for i in 1 2 3 4 5 6
do
    sample=${samples[${i}]}
    echo ${sample}
    #Map the reads
    /usr/local/tophat-2.0.8/tophat -p 20 -G ${annotation} -o /home/shahrokh/bio720/FinalProject/mapped/${sample} ${reference} \
     /home/shahrokh/bio720/FinalProject/trimmed/${sample}_1P.fastq.gz \
     /home/shahrokh/bio720/FinalProject/trimmed/${sample}_2P.fastq.gz \
    #sort out the mapped reads based on names
    samtools sort -n \
    ~/bio720/FinalProject/mapped/${sample}/accepted_hits.bam ~/bio720/FinalProject/mapped/${sample}/accepted_hits_sorted.bam
	
    #Count the number of reads mapping to each feature using HTSeq
    htseq-count --format bam --stranded no --type gene --idattr ID --order name --mode intersection-nonempty \
    ~/bio720/FinalProject/mapped/${sample}/accepted_hits_sorted.bam ${annotation} > ~/bio720/FinalProject/counts/${sample}_htseq_counts.txt
done
