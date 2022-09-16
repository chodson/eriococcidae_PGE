First step: checked read quality with FastQC

## Transcriptome assemblies

### Fastp to trim reads
example command for one library below
```
fastp -i 0_reads/3538F4.1.r1.fastq.gz -I 0_reads/3538F4.1.r2.fastq.gz -o 2_trim/3538F4.1.r1.trim.fastq.gz -O 2_trim/3538F4.1.r2.trim.fastq.gz --cut_by_quality5 --cut_by_quality3 --cut_window_size 4 --cut_mean_quality 20 --trim_poly_g --html 2_trim/3538F4.1.html --thread 4
```
- read quality looked good after trimming


### Trinity to assemble transcriptome
<https://github.com/trinityrnaseq/trinityrnaseq/wiki>
below is an example of the command used for C. campanidorsalis, the same parameters were used for C. echiniformis libraries. The sample file for C. echiniformis included all samples from the three families

```
/ceph/software/trinity/trinityrnaseq-Trinity-v2.8.4/Trinity --seqType fq --max_memory 100G --samples_file 3_trinity/ccamp/sample_file.txt --SS_lib_type RF --CPU 30  --output 3_trinity/ccamp/trinity --workdir /scratch/chodson/trinity --full_cleanup
```

sample file- includes trimmed reads
1f.ccamp.3538F4	ccamp.3538F4.f.1	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.1.r1.trim.fastq.gz	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.1.r2.trim.fastq.gz
2ccamp.3538F4	ccamp.3538F4.f.2	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.2.r1.trim.fastq.gz	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.2.r2.trim.fastq.gz
3ccamp.3538F4	ccamp.3538F4.f.3	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.3.r1.trim.fastq.gz	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.3.r2.trom.fastq.gz
m1ccamp.3538F4	ccamp.3538F4.m1	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.m1.r1.trim.fastq.gz	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.m1.r2.trim.fastq.gz
cm2camp.3538F4	ccamp.3538F4.m2	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.m2.r1.trim.fastq.gz	/data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.m2.r2.trim.fastq.gz	


### Use RSEM to get FPKM values for each isoform and filter out isoforms with very little/ no expression 
<https://github.com/deweylab/RSEM>
- other papers filter in a similar way (Bast 2018, they use RKPM) to exclude isoforms which don't have evidence they are expressed (because Trinity predicts too many isoforms)
- I did this for each lib separately and kept any transcripts that had FPKM 0.5 or higher in at least one lib.
- Used FPKM cutoff of 0.5 since the goal was to filter out transcripts with no or very little expression
- used seqtk and filtered transcriptome by IDs that met criteria


```
qsub -o logs -e logs -cwd -N rsem -V -pe smp64 1 -b yes 'rsem-prepare-reference --bowtie2 transcriptomes/ccamp.transcriptome.fasta transcriptomes/ccamp.transcriptome'
 
qsub -o logs -e logs -cwd -N rsem -V -pe smp64 1 -b yes 'rsem-prepare-reference --bowtie2 transcriptomes/cech.transcriptome.fasta transcriptomes/cech.transcriptome'
```
example of the command used for one C. camp and one C. ech library. Done for all samples
```
qsub -o logs -e logs -cwd -N rsem -V -pe smp64 16 -b yes 'rsem-calculate-expression --paired-end --bowtie2 -p 16 /data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.1.r1.trim.fastq.gz /data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3538F4.1.r2.trim.fastq.gz transcriptomes/ccamp.transcriptome 4_rsem/ccamp3538F4.f1'
```
```
qsub -o logs -e logs -cwd -N rseme1 -V -pe smp64 16 -b yes 'rsem-calculate-expression --paired-end --bowtie2 -p 16 /data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3571F5.r1.trim.fastq.gz /data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3571F5.r2.trim.fastq.gz transcriptomes/cech.transcriptome 4_rsem/cech3571F5'
```

#### Filtering out isoforms with FPKM>0.5
example of the command used for one C. camp and one C. ech library. Done for all samples
```
awk '$7 > 0.5' 4_rsem/ccamp3538F4.f1.isoforms.results > 4_rsem/ccamp3538F4.f1.isoforms.filtered2.results 

awk '$7 > 0.5' 4_rsem/cech3571F5.isoforms.results > 4_rsem/cech3571F5.isoforms.filtered2.results 
```

```
cat ccamp*filtered2.results | cut -f1 | sort | uniq > c.camp2.filteredFPKM.id.txt
cat cech*filtered2.results | cut -f1 | sort | uniq > c.ech2.filteredFPKM.id.txt
```
#### Filtering out transcripts with FPKM<0.5 from transcriptome 
```
qsub -o logs -e logs -cwd -N rseme3 -V -pe smp64 16 -b yes 'seqtk subseq transcriptomes/ccamp.transcriptome.fasta 4_rsem/c.camp2.filteredFPKM.id.txt > transcriptomes/ccamp.transcriptome.filtered.fasta'
qsub -o logs -e logs -cwd -N seqtk -V -pe smp64 1 -b yes 'seqtk subseq transcriptomes/cech.transcriptome.fasta 4_rsem/c.ech2.filteredFPKM.id.txt > transcriptomes/cech.transcriptome.filtered.fasta'
```

### Used trinity script get_longest_isoform_seq_per_trinity_gene.pl to filter out longest orf from remaining transcripts

#### Getting longest orf of remaining transcripts
```
qsub -o logs -e logs -cwd -N longiso -V -pe smp64 1 -b yes '/ceph/users/kjaron/.conda/envs/javabased/./opt/trinity-2.8.5/util/misc/get_longest_isoform_seq_per_trinity_gene.pl transcriptomes/ccamp.transcriptome.filtered.fasta > transcriptomes/ccamp.transcriptome.filtered.longiso.fasta'
qsub -o logs -e logs -cwd -N longiso -V -pe smp64 1 -b yes '/ceph/users/kjaron/.conda/envs/javabased/./opt/trinity-2.8.5/util/misc/get_longest_isoform_seq_per_trinity_gene.pl transcriptomes/cech.transcriptome.filtered.fasta > transcriptomes/cech.transcriptome.filtered.longiso.fasta'
```
```
qsub -o logs -e logs -cwd -N td.lorf -V -pe smp64 1 -b yes 'TransDecoder.LongOrfs -t transcriptomes/ccamp.transcriptome.filtered.longiso.fasta --output_dir 5_transdecoder2/ccamp/'
qsub -o logs -e logs -cwd -N td.lorf -V -pe smp64 1 -b yes 'TransDecoder.LongOrfs -t transcriptomes/cech.transcriptome.filtered.longiso.fasta --output_dir 5_transdecoder2/cech/'
```

### Using transdecoder to retain ORFs with homology to something.
- used transdecoder longorf with default settings
- blastp with swissprot database and parameters: -max_target_seqs 1 -outfmt 6 -evalue 1e-5
- hmmscan with default parameters and database /ceph/software/blast_db/pfam-31.0/Pfam-A.hmm
- transdecoder predict: --single_best_only parameter and hmmscan and blastp outputs.

```
qsub -o logs -e logs -cwd -N blast.td -V -pe smp64 16 -b yes 'blastp -query 5_transdecoder2/ccamp/longest_orfs.pep  -db /ceph/software/blast_db/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 > 5_transdecoder2/ccamp/ccamp.vs.uniprot_sprot.blastp.outfmt6'
qsub -o logs -e logs -cwd -N hmm.td -V -pe smp64 16 -b yes 'hmmscan --cpu 16 --domtblout 5_transdecoder2/ccamp/ccamp.vs.pfam.domtblout /ceph/software/blast_db/pfam-31.0/Pfam-A.hmm 5_transdecoder2/ccamp/longest_orfs.pep'
qsub -o logs -e logs -cwd -N blast.td -V -pe smp64 16 -b yes 'blastp -query 5_transdecoder2/cech/longest_orfs.pep  -db /ceph/software/blast_db/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16 > 5_transdecoder2/cech/cech.vs.uniprot_sprot.blastp.outfmt6'
qsub -o logs -e logs -cwd -N hmm.td -V -pe smp64 16 -b yes 'hmmscan --cpu 16 --domtblout 5_transdecoder2/cech/cech.vs.pfam.domtblout /ceph/software/blast_db/pfam-31.0/Pfam-A.hmm 5_transdecoder2/cech/longest_orfs.pep'
```
```
qsub -o logs -e logs -cwd -N pred.td -V -pe smp64 1 -b yes 'TransDecoder.Predict --single_best_only -t transcriptomes/ccamp.transcriptome.fasta --retain_pfam_hits 5_transdecoder2/ccamp/ccamp.vs.pfam.domtblout --retain_blastp_hits 5_transdecoder2/ccamp/ccamp.vs.uniprot_sprot.blastp.outfmt6 --output_dir 5_transdecoder2/ccamp/'
qsub -o logs -e logs -cwd -N pred.td -V -pe smp64 1 -b yes 'TransDecoder.Predict --single_best_only -t transcriptomes/cech.transcriptome.fasta --retain_pfam_hits 5_transdecoder2/cech/cech.vs.pfam.domtblout --retain_blastp_hits 5_transdecoder2/cech/cech.vs.uniprot_sprot.blastp.outfmt6 --output_dir 5_transdecoder2/cech/'
```

## SNP calling

### Mapping reads to transcriptome with Bowtie2
<http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>

#### Making bowtie2 index with filtered transdecoder output
```
qsub -cwd -e logs -o logs -N bowtie_index -V -pe smp64 1 -b yes 'bowtie2-build transcriptomes/ccamp.transcriptome.fasta.transdecoder.filtered.longiso.cds transcriptomes/ccamp.transcriptome.fasta.transdecoder.filtered.longiso'
qsub -cwd -e logs -o logs -N bowtie_index -V -pe smp64 1 -b yes 'bowtie2-build transcriptomes/cech.transcriptome.fasta.transdecoder.filtered.longiso.cds transcriptomes/cech.transcriptome.fasta.transdecoder.filtered.longiso'
```

#### Mapping reads

- example for 3571F5.m4 sample
```
qsub -cwd -e logs -o logs -N bowtie2 -V -pe smp64 16 -b yes 'bowtie2 --rg-id cech3571F5M4_cech3571F5M4_lane1 --rg LB:cech3571F5M4_cech3571F5M4 --rg PL:ILLUMINA --rg PU:cech3571F5M4_cech3571F5M4_lane1 --rg SM:cech3571F5M4 -p 16 -1 /data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3571F5.m4.r1.trim.fastq.gz -2 /data/ross/mealybugs/analyses/Cystococcus.rna/2_trim/3571F5.m4.r2.trim.fastq.gz -x transcriptomes/cech.transcriptome.fasta.transdecoder.filtered.longiso | samtools sort -T /scratch/chodson/bowtie -@8 -O bam - > /scratch/chodson/3571F5.m4.vs.filteredlongiso.bowtie.ccamp.bam && rsync --remove-source-files /scratch/chodson/3571F5.m4.vs.filteredlongiso.bowtie.ccamp.bam 6_mapping'
```
- merging F mapping from different lanes into one bam
```
qsub -cwd -e logs -o logs -N st.sort -V -pe smp64 8 -b yes 'samtools merge -@ 8 6_mapping/3538F4.vs.filteredlongiso.bowtie.ccamp.bam 6_mapping/3538F4.1.vs.filteredlongiso.bowtie.ccamp.bam 6_mapping/3538F4.2.vs.filteredlongiso.bowtie.ccamp.bam 6_mapping/3538F4.3.vs.filteredlongiso.bowtie.ccamp.bam'
```
#### Marking duplicates in bam files (for all bam files)
```
qsub -cwd -e logs -o logs -N st.sort -V -pe smp64 8 -b yes 'samtools sort -n -@ 8 -T /scratch/chodson/ -o 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.namesort.bam 6_mapping/3538F4.vs.filteredlongiso.bowtie.ccamp.bam && samtools fixmate -m 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.namesort.bam 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.fixmate.bam && samtools sort -@ 8 -T /scratch/chodson/ -o 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.positionsort.bam 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.fixmate.bam && samtools markdup -T /scratch/chodson/ 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.positionsort.bam 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.markdup.bam'
```

#### Freebayes to call variants

```
qsub -cwd -e logs -o logs -N freebayes -V -pe smp64 1 -b yes 'freebayes -f transcriptomes/ccamp.transcriptome.fasta.transdecoder.filtered.longiso.cds -b 7_markdups/3538F4.m1.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3538F4.m2.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.markdup.bam --standard-filters --vcf /scratch/chodson/ccamp.filtered.longiso.variants.vcf -C 5 --min-coverage 10 -w && rsync --remove-source-files /scratch/chodson/ccamp.filtered.longiso.variants.vcf 8_freebayes/'
```
```
qsub -cwd -e logs -o logs -N freebayes -V -pe smp64 1 -b yes 'freebayes -f transcriptomes/cech.transcriptome.fasta.transdecoder.filtered.longiso.cds -b 7_markdups/3572F4.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3572F4.m3.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3572F4.m2.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3571F6.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3571F6.m.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3571F6.m2.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3571F5.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3571F5.m4.vs.filteredlongiso.bowtie.ccamp.markdup.bam 7_markdups/3571F5.m3.vs.filteredlongiso.bowtie.ccamp.markdup.bam --standard-filters --vcf /scratch/chodson/cech.filtered.longiso.variants.vcf -C 5 --min-coverage 10 -w && rsync --remove-source-files /scratch/chodson/cech.filtered.longiso.variants.vcf 8_freebayes/'
```


#### Filtering variants
Filtered in steps because I was checking whether what I was getting after each step make sense. In the end...
Filtering criteria:
- Only retained SNPs 
- with a depth greater than 10
- and a Quality score greater than 20
- and an allelic balance AB between 0.2-0.8

Commands:

```
qsub -cwd -e logs -o logs -N snpfilt -V -pe smp64 1 -b yes 'vcffilter -f "LEN = 1 & TYPE = snp & NUMALT = 1" 8_freebayes/cech.filtered.longiso.variants.vcf > 9_filtersnp/cech.filtered.longiso.variants.snp.vcf && vcfstats 9_filtersnp/cech.filtered.longiso.variants.snp.vcf > snp.stats.cech.snponly'
```
```
qsub -cwd -e logs -o logs -N snpfilt -V -pe smp64 1 -b yes 'vcffilter -f "DP > 10 & SAF > 2 & SAR > 2" 9_filtersnp/cech.filtered.longiso.variants.snp.vcf > 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.vcf && vcfstats 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.vcf > snp.cech.stats.snponly.dp10'
```
```
qsub -cwd -e logs -o logs -N snpfilt -V -pe smp64 1 -b yes 'vcffilter -f "QUAL > 20" 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.vcf > 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.vcf && vcfstats 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.vcf > snp.cech.stats.snponly.dp10.q && vcffilter -f "QUAL > 20" 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.vcf > 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.vcf && vcfstats 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.vcf > snp.camp.stats.snponly.dp10.q'
```
```
qsub -cwd -e logs -o logs -N snpfilt -V -pe smp64 1 -b yes 'vcffilter -f "AB > 0.2" 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.vcf > 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.ab.vcf && vcfstats 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.ab.vcf > snp.cech.stats.snponly.dp10.q.ab && vcffilter -f "AB > 0.2" 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.vcf > 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab.vcf && vcfstats 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab.vcf > snp.ccamp.stats.snponly.dp10.q.ab'
```
```
qsub -cwd -e logs -o logs -N snpfilt -V -pe smp64 1 -b yes 'vcffilter -f "AB < 0.8" 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.ab.vcf > 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.ab2.vcf && vcfstats 9_filtersnp/cech.filtered.longiso.variants.snp.dp10.q.ab2.vcf > snp.cech.stats.snponly.dp10.q.ab2 && vcffilter -f "AB < 0.8" 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab.vcf > 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab2.vcf && vcfstats 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab2.vcf > snp.ccamp.stats.snponly.dp10.q.ab2'
```
#### ASE readcounter commands
- below is the commands for C. campanidorsalis. The same criteria were used for C. echiniformis
```
qsub -cwd -e logs -o logs -N ase -V -pe smp64 1 -b yes 'gatk ASEReadCounter -R transcriptomes/ccamp.transcriptome.fasta.transdecoder.filtered.longiso.fasta -I 7_markdups/3538F4.m1.vs.filteredlongiso.bowtie.ccamp.markdup.bam -O 11_counting_reads/transcriptome.vs.3538.m1.csv -V 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab2.vcf --min-depth-of-non-filtered-base 20 --min-base-quality 20 --min-mapping-quality 40'
qsub -cwd -e logs -o logs -N ase -V -pe smp64 1 -b yes 'gatk ASEReadCounter -R transcriptomes/ccamp.transcriptome.fasta.transdecoder.filtered.longiso.fasta -I 7_markdups/3538F4.m2.vs.filteredlongiso.bowtie.ccamp.markdup.bam -O 11_counting_reads/transcriptome.vs.3538.m2.csv -V 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab2.vcf --min-depth-of-non-filtered-base 20 --min-base-quality 20 --min-mapping-quality 40'
qsub -cwd -e logs -o logs -N ase -V -pe smp64 1 -b yes 'gatk ASEReadCounter -R transcriptomes/ccamp.transcriptome.fasta.transdecoder.filtered.longiso.fasta -I 7_markdups/3538F4.vs.filteredlongiso.bowtie.ccamp.markdup.bam -O 11_counting_reads/transcriptome.vs.3538.f.csv -V 9_filtersnp/ccamp.filtered.longiso.variants.snp.dp10.q.ab2.vcf --min-depth-of-non-filtered-base 20 --min-base-quality 20 --min-mapping-quality 40'
```

#### next steps

- Computed for each sample the level of heterozygous allele expression.
- I did this in R with the script:
`~/projects/cystococcus_snp/cystococcus_other/bits_for_paper.R`