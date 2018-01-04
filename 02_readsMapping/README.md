# 二代测序 reads 的 mapping
二代测序 reads 的比对分两步：
  1. 用 bwa 将 reads 比对到参考基因组
  2. 用 picard 去除 PCR 造成的重复 reads
 
最终得到比对后的 CRAM 文件，这是一种BAM的压缩格式，在 samtools 给的[基准测试](http://www.htslib.org/benchmarks/CRAM.html)中，CRAM 大小约为 BAM 的一半。
# 01 bwa比对
[bwa](https://github.com/lh3/bwa)软件实现了三种比对算法 BWA-backtrack, BWA-SW 和 BWA-MEM。第一种算法适用于长度在 100bp 以下的 reads，后两种算法适用于70bp至数M的长 reads。BWA-MEM 是最新的算法，70bp以上的 reads 用 BWA-MEM 就好，70b p以下的用 BWA-backtrack。
## 01.1 bwa-mem
```
bwa \
    mem \
    -t 4 \
    -R '@RG\tID:sample1_library1_lane1\tLB:library1\tPL:ILLUMINA\tSM:sample1' \
    reference.fa \
    sample1_library1_lane1_1.clean.fq.gz \
    sample1_library1_lane1_2.clean.fq.gz | \
samtools \
    fixmate \
    -r \
    - - | \
samtools \
    sort \
    --reference reference.fa \   
    -o sample1_library1_lane1.bam \
    --output-fmt BAM \
    -
```
这一步对reads进行了比对(mem)，并去掉未比对到参考基因组上的reads(fixmate -r)，最后又把BAM文件根据参考基因组的坐标位置进行了排序(sort)。
- 注:
1. -t 为使用的线程数，-R为BAM文件的header，需要根据实际样本信息去修改。
2. 通常来说，这步得到的文件还需要使用 picard 去除 PCR 重复，如果后续不使用 picard 去重（或者只使用 samtools 去重） 的话，这一步可以直接生成 cram 文件，即把 --output-fmt BAM 改为 CRAM ，同时把 -o 的文件名后缀改为.cram
3. picard去重复不需要对BAM文件建索引，所以如果后续还要去重的话，就不需要对这步的结果建索引。
4.这里使用了samtools fixmate来修复一些bwa直接输出的FLAG（以更好地兼容下游软件），同时加了-r选项来去除未比对上的reads以减少文件大小，视情况可以不加
## 01.2 bwa-backtrack
如果使用bwa backtrack比对需要两步，第一步的命令是bwa aln，第二步有两种，单端测序用bwa samse，双端测序用bwa sampe
### 01.2.1 bwa aln
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
    -o sample1_library1_lane1.bam \
    --output-fmt BAM \
    -
```
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
    -o sample1_library1_lane1.bam \
    --output-fmt BAM \
    -
```

## 02 picard去重
```
java -Djava.io.tmpdir=<tmpdir> \
    -Xmx10g \
    -jar picard.jar \
    MarkDuplicates \
    I=sample1_library1_lane1.bam \
    I=sample1_library1_lane2.bam \
    O=/dev/stdout \
    M=sample1.marked_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT \
    COMPRESSION_LEVEL=0 | \
samtools \
    view \
    -T reference.fa \
    -C \
    -o sample1.cram \
    -

samtools \
    index \
    sample1.cram
```
这一步是把同一个体的多个BAM文件同时去重，最终每个个体生成一个CRAM文件，最后建立索引。
- 注：
1. -Djava.io.tmpdir=<tmpdir>指定了输出临时文件的目录，使用时把<tmpdir>替换为自己有写入权限的目录
2. COMPRESSION_LEVEL=0是为了便于samtools直接读取并转为cram格式
3. REMOVE_DUPLICATES=true是直接从BAM文件里删掉了重复的reads（而不是只做标记）
4. 使用脚本produce_picard_dedup.py可根据BAM文件的header批量生成去重的shell脚本。
