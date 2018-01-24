# bioNotes

### one liner

fastq total base\n
zcat 1.fq.gz | awk '(NR % 4) == 2 {a+=length($1)}; END{print a}'
