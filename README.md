# bioNotes

### one liner

- fastq total base
zcat 1.fq.gz | awk '(NR % 4) == 2 {a+=length($1)}; END{print a}'
