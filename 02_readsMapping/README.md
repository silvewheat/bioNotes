# 二代测序reads的mapping
[bwa](https://github.com/lh3/bwa)软件实现了三种比对算法BWA-backtrack, BWA-SW 和 BWA-MEM。第一种算法适用于长度在100bp以下的reads，后两种算法适用于70bp至数M的长reads。BWA-MEM是最新的算法，70bp以上的reads用BWA-MEM就好，70bp以下的用BWA-backtrack。
## bwa-backtrack
bwa backtrack比对需要两步，第一步的命令是bwa aln，第二步有两种，单端测序用bwa samse，双端测序用bwa sampe
### bwa aln
```
    bwa \
        aln \
        -t 8 \
        -l 1024 \
        reference.fa \
        sample1_library1_lane1_1.clean.fq.gz \
        > sample1_library1_lane1_1.sai
```
其中-t指定线程数，-l 1024意味着禁用种子。
单端测序的reads只生成一个.sai文件，双端测序_1和_2的reads分别生成对应的.sai文件。
produce_bwa_aln.py可用于批量生成bwa aln的shell脚本。
### bwa samse
单端测序的reads经过bwa aln得到的.sai文件作为bwa samse的输出文件，生成最终的比对结果文件。
结合管道符可以一步生成排序后的结果。
```
bwa \
    samse \
    -r '@RG\tID:sample1_library1_lane1\tLB:library1\tPL:ILLUMINA\tSM:sample1' \
    reference.fa \
    sample1_library1_lane1_1.sai \
    sample1_library1_lane1_1.clean.fq.gz | \
samtools \
    fixmate \
    -r \
    - - | \
samtools \
    sort \
    --reference reference.fa \   
    -o sample1_library1_lane1.cram \
    --output-fmt CRAM \
    -
```
这里使用了samtools fixmate来修复一些bwa直接输出的FLAG（以兼容下游软件），同时加了-r选项来去除未比对上的reads以减少文件大小。
最终输出的格式选择了CRAM文件，这是一种BAM的压缩格式，在samtools给的[基准测试](http://www.htslib.org/benchmarks/CRAM.html)中，CRAM大小约为BAM的一半。
produce_bwa_samse_sampe.py可用于批量生成bwa samse和bwa sampe的shell脚本。
### bwa sampe
与bwa sampe类似
```
bwa \
    sampe \
    -r '@RG\tID:sample1_library1_lane1\tLB:library1\tPL:ILLUMINA\tSM:sample1' \
    reference.fa \
    sample1_library1_lane1_1.sai \
    sample1_library1_lane1_2.sai \
    sample1_library1_lane1_1.clean.fq.gz \
    sample1_library1_lane1_2.clean.fq.gz | \
samtools \
    fixmate \
    -r \
    - - | \
samtools \
    sort \
    --reference reference.fa \   
    -o sample1_library1_lane1.cram \
    --output-fmt CRAM \
    -
```
