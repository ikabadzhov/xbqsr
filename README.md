## SYCL-GAL
#### In-memory, single-Threaded implementation of GATK best practice pipeline

##### Note:

The current version support the following stages:
+ sorting (by coordinates only)
+ marking duplicates (with optical duplicates count)
+ bqsr 
  + base ricalibrator
  + apply bqsr

Supported I/0: 
+ sam 
+ bam (not fully tested)

### Example usage

```
./GenomicPipeline --input aligned_reads.sam --output output.sam --sorting-order coordinates --mark-duplicates --metric-file metrics.txt --bqsr recall.txt --reference Homo_sapiens_assembly38.fasta --known-sites Homo_sapiens_assembly38.known_indels.vcf
```

