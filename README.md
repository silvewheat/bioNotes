# bioNotes

### one liner

fastq total base :

    zcat 1.fq.gz | awk '(NR % 4) == 2 {a+=length($1)}; END{print a}'

extract specific reads from fastq file according to reads name :

    zcat a.fastq.gz | awk 'BEGIN{RS="@";FS="\n"}; $1~/readsName/{print $2; exit}'

count missing sample in vcf file per line:

    bcftools query -f '[%GT\t]\n' a.bcf |  awk '{miss=0};{for (x=1; x<=NF; x++) if ($x=="./.") {miss+=1}};{print miss}' > nmiss.count
