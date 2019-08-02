## Pre-processing

Preliminary analysis of soil samples were done using [FAST](https://github.com/ZeweiSong/FAST), reads were trimmed and processed using cutadapt and vsearch.
Before processing with **FAST**, files were extracted from gz format since there was an error in one of the scripts that failed to recognized the encoding.  The first step is to generate a quick map based qiime guidelines.

```
#Mapping file for forward reads
python ~/bin/FAST/fast.py -generate_mapping -i reads_1 -o read1_map.txt
#Mapping file for reverse reads
python ~/bin/FAST/fast.py -generate_mapping -i reads_2 -o read2_map.txt
```

The mapping files should follow all the parameters designated by [qiime developers](http://qiime.org/scripts/add_qiime_labels.html).  The next step is to add labels to the fastq files to match the mapping file. `-t 4` _this flag indicates a parallel processing using 4 threads_.

```
#Add labels to fastq files - forward reads
python ~/bin/FAST/fast.py -add_labels -m read1_map.txt -i reads_1 -o read1_labeled -t 8
#Add labels to fastq files - reverse reads
python ~/bin/FAST/fast.py -add_labels -m read2_map.txt -i reads_2 -o read2_labeled -t 8
```

After labeling reads per file using the sample name, which refers to the mapping file, reads could be merged into a single fastq file for downstream analyses.  This eases the downstream anlayses reducing the number of commands per file.

Merge all labeled sequences:

```
#Merge forward reads into a single fastq file
python ~/bin/FAST/fast.py -merge_seqs -i read1_labeled -o read1.fastq
#Merge reverse reads into a single fastq file
python ~/bin/FAST/fast.py -merge_seqs -i read2_labeled -o read2.fastq
```

To further clean the reads, using cutadapt, illumina adapters (or whatever adapters were used) could be removed from reads and at the same time remove any reads that are lest than 50 bp.  This is important before assembling the reads into a single sequence.

In our case, these are the Illumina overhang adapterused:

- Forward overhang 5’ GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG-[locus-specific sequence]
- Reverse overhang 5’ GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-[locus-specific sequence]

```
cutadapt -a GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG -A GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
          -o read1.cut.fastq \
          -p read2.cut.fastq read1.fastq read2.fastq \
          -m 50 -j 6
```
These are the results after removing the adapters.  Despite that the sequencing center removed adpaters, there was an small percentage still left in the reads.

```
=== Summary ===
Total read pairs processed:          2,821,426
  Read 1 with adapter:                  17,217 (0.6%)
  Read 2 with adapter:                   8,922 (0.3%)
Pairs that were too short:                   1 (0.0%)
Pairs written (passing filters):     2,821,425 (100.0%)

Total basepairs processed: 1,698,498,452 bp
  Read 1:   849,249,226 bp
  Read 2:   849,249,226 bp
Total written (filtered):  1,698,413,769 bp (100.0%)
  Read 1:   849,194,042 bp
  Read 2:   849,219,727 bp
```

To further clean the reads and focus only on the variable region of ITS, primers ITS1f and ITS could be reomved:

```
cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCTGCGTTCTTCATCGATGC \
         -o read1.test.fastq -p read2.test.fastq \
         read1.fastq read2.fastq -e 0.1 \
         --match-read-wildcards -m 50 -j 10
```

This results in:
```
#Result
=== Summary ===

Total read pairs processed:          2,821,426
  Read 1 with adapter:               2,796,053 (99.1%)
  Read 2 with adapter:                       0 (0.0%)
Pairs that were too short:                  11 (0.0%)
Pairs written (passing filters):     2,821,415 (100.0%)
```
After removing the Illumina adapters, the reads were assembled using [PEAR]():

```
#The flag '-j 4' dictates the number threads used for processing the reads
pear -f read1.cut.fastq -r read2.cut.fastq -o merge.pear -j 10
```

The results for merging the reads are below, most of the reads were assembled, keeping 97.6% of thre data.  Only a small percentage did not assemble (1.139%), these unassemble reads have reduced quality or the reads do not pair.

```
Assembled reads ...................: 2,784,935 / 2,821,415 (98.707%)
Discarded reads ...................: 0 / 2,821,415 (0.000%)
Not assembled reads ...............: 36,480 / 2,821,415 (1.293%)
Assembled reads file...............: merge.pear.assembled.fastq
Discarded reads file...............: merge.pear.discarded.fastq
Unassembled forward reads file.....: merge.pear.unassembled.forward.fastq
Unassembled reverse reads file.....: merge.pear.unassembled.reverse.fastq
```

Using vsearch, low quality sequences can be filtered, using a max expected error (<1):
- --fastq_maxee REAL  maximum expected error value for filter
- --fasta_width INT   width of FASTA seq lines, 0 for no wrap (80)

```
vsearch --fastq_filter merge.pear.assembled.fastq \
        --fastq_maxee 1 \
        --fastaout merge.pear.maxee1.fasta \
        --fasta_width 0 --threads 10
```
It resulted in 0.6% of the assembled reads discarded due to low quality.

```
2767437 sequences kept (of which 0 truncated), 17498 sequences discarded.
```

#De-replication and OTU clustering

Using vsearch through FAST, reads can be dereplicated to ease the downstream analyses.

```
python ~/bin/FAST/fast.py -dereplicate -i merge.pear.maxee1.fasta -o raw.qc.derep -t 10
```
A file `raw.qc.derep.txt` containing the information for the OTU map was written. This file is basedon clusters of sequences with  100% similarity.  A second file `raw.qc.derep.fasta` is also created that contains the unique sequences after this dereplication step. The summary of the results is presented below:

```
Original OTU map:
	 OTU=2145966 (Total Sequences=2767437, Max=2845, Min=1, Ave=1)
Filtered OTU map:
	 OTU=179850 (Total Sequences=801321, Max=2845, Min=2, Ave=4)
```
Singletons can be removed at this step, since those are artifacts or they can be removed later based on their presence across different samples.  In this case, I will remove singletons here:

```
python ~/bin/FAST/fast.py -filter_otu_map \
                                  -i raw.qc.test.derep.txt \
                                  -o raw.qc.derep.size2.txt -min_size 2
```
This is the resulting output:

```
Original OTU map:
	 OTU=2145966 (Total Sequences=2767437, Max=2845, Min=1, Ave=1)
Filtered OTU map:
	 OTU=179850 (Total Sequences=801321, Max=2845, Min=2, Ave=4)
```

New sequence file after removing singletons, and adding size label:
```
python ~/bin/FAST/fast.py -pick_seqs -i raw.qc.test.derep.fasta \
                                  -map raw.qc.derep.size2.txt \
                                  -o raw.qc.derep.size2.fasta -sizeout
```
Here the `-sizeout` option will add a size annotation to each sequence `(;size=XXX)`. This is require by VSEARCH for chimera checking and OTU clustering.

One last thing is to check for chimeras before an OTU analysis is done. 

```
vsearch --uchime_ref raw.qc.derep.size2.fasta \
        --nonchimeras raw.qc.derep.size2.uchime.fasta \
        --db /data2/unite_qiime/sh_refs_qiime_ver8_99_s_02.02.2019.fasta \
        --sizeout --fasta_width 0 --thread 10
```

The chimera content is low, around 6.9% taking into account abundance information.

```
Found 13610 (7.6%) chimeras, 157581 (87.6%) non-chimeras,
and 8659 (4.8%) borderline sequences in 179850 unique sequences.
Taking abundance information into account, this corresponds to
55641 (6.9%) chimeras, 708496 (88.4%) non-chimeras,
and 37184 (4.6%) borderline sequences in 801321 total sequences.
```

Finally, the OTU clustering could be done setting up a 97% similarity using flag `-id`.

```
vsearch --cluster_size raw.qc.derep.size2.uchime.fasta \
        --centroids raw.qc.vsearch.fasta \
        --fasta_width 0 -id 0.97 \
        --sizein --uc raw.qc.uc.txt --threads 10
```

There are 943 OTUs, and there are no singletons, since those were removed previously.
```
Clusters: 943 Size min 2, max 117604, avg 167.1
Singletons: 0, 0.0% of seqs, 0.0% of clusters
```

Generating an OTU amp using UC file `raw.qc.uc.txt`. The OTU sequences are in the `raw.qc.vsearch.fasta` file.

```
python ~/bin/FAST/fast.py -parse_uc_cluster -i raw.qc.uc.txt -o raw.qc.vsearch.txt
```

In FAST pipeline OTU map and sequences are combined into a single file. It is in JSON format that can be readily read into a Python dictionary. 

Combine the Dereplicate map and sequences:
```
python ~/bin/FAST/fast.py -generate_fast_map \
                                  -map raw.qc.derep.size2.txt \
                                  -seq raw.qc.derep.size2.uchime.fasta \
                                  -o fast.derep.txt -derep
```

Combine the OTU map and sequences:
```
python ~/bin/FAST/fast.py -generate_fast_map \
                                  -map raw.qc.vsearch.txt \
                                  -seq raw.qc.vsearch.fasta \
                                  -o fast.otu.txt -otu
```

Combine to FAST derep map and OTU map into a single hybrid:
```
python ~/bin/FAST/fast.py -combine_fast_map \
                                  -derep_map fast.derep.txt \
                                  -otu_map fast.otu.txt -o fast.hybrid.txt
```

Rename the OTUs, so them will start with `OTU_`:
```
python ~/bin/FAST/fast.py -rename_otu_map -fast_map fast.hybrid.txt -o fast.hybrid.otu.txt
```

Generate the OTU table from the FAST hybrid map, along with the representative sequences:
```
python ~/bin/FAST/fast.py -make_otu_table \
                                  -fast_map fast.hybrid.otu.txt \
                                  -o otu_table.txt -rep rep_seq.fasta
```

# Assigning taxonomy

We can now align our OTU sequences to the UNITE database using Blast. 

```
parallel_assign_taxonomy_blast.py -O 16 -t /data2/unite_qiime/sh_taxonomy_qiime_ver8_dynamic_s_02.02.2019.txt \
                                                                  -r /data2/unite_qiime/sh_refs_qiime_ver8_dynamic_s_02.02.2019.fasta \
                                                                  -i rep_seq.fasta \
                                                                  -o blast_all_assigned_taxonomy
```

# Make otu table biom file
Using qiime1, a biom file is created and the taxonomy and sample data could be added:
````
#Create biom file
biom convert --table-type="OTU table" -i otu_table_mapping.txt -o otu_table_mapping.biom --to-json

#adding taxonomy
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy \
				  --observation-metadata-fp blast_all_assigned_taxonomy/rep_seq_tax_assignments.txt \
				  -i otu_table_TomAR.biom -o otu_table_TomAR.tax.biom --output-as-json
				  
#adding sample data
biom add-metadata -m read_map.txt -i otu_table_TomAR.tax.biom -o otu_table_TomAR.tax.sample.biom --output-as-json
```

# Filter NTC from samples
```
#create Biom for negative controls (NTCs)
filter_samples_from_otu_table.py -i otu_table_TomAR.tax.sample.biom -o otu_table_NTC.biom -m read_map.txt -s "Samples:NTC"

#Biom without NTC samples
filter_samples_from_otu_table.py -i otu_table_TomAR.tax.sample.biom -o otu_table_TomAR1.tax.sample.biom -m read_map.txt -s "Samples:AR.Tom"

#remove non zeros from NTCs
filter_otus_from_otu_table.py -i otu_table_NTC.biom -o otu_table_NTC_Non.biom -n 1

#Convert NTC_biom file to a tab-separated txt table
biom convert -i otu_table_NTC_Non.biom -o otus_to_remove.txt --to-tsv

#Filter NTC OTUS using txt file from NTC samples
filter_otus_from_otu_table.py -i otu_table_TomAR1.tax.sample.biom -o otu_table_TomasAR.clean.biom -e otus_to_remove.txt

#Filter non-zero samples
filter_otus_from_otu_table.py -i otu_table_TomasAR.clean.biom -o otu_table_TomasAR.clean1.biom -n 1

#Check counts in final biom file
biom summarize-table -i otu_table_TomasAR.clean1.biom
```

