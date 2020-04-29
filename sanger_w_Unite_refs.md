# Including sanger sequences in the UNITE database

The most recent version of the Unite (Unite SH v8.2) was downloaded from [unite.ut.ee](https://unite.ut.ee/repository.php) and used to classify sanger sequences using blast or RDP classifier

```
wget -O Unite_fungi.tar.gz https://files.plutof.ut.ee/public/orig/1E/66/1E662B6EB320312A61E7E3218327F34C7DB09CFF8E4686A89EF47886822DA6AB.gz

gunzip Unite_fungi.tar.gz

tar -xvf Unite_fungi.tar
```

The fasta file with the dynamic assignment of the SH was used to generate a blast database

```
makeblastdb -in sh_general_release_dynamic_s_04.02.2020.fasta -dbtype 'nucl' -out unite_040220
```

Then the sanger sequences were blasted against the blast database from the previous step

```
blastn -db unite_040220 -query Milani_Sanger_GoodSequences.fasta -out tax_Tomas_sanger.txt -outfmt "6 std stitle"
```

Now, the output file will be parsed to generate the taxonomy file
```
cat tax_Tomas_sanger.txt | cut -d$'\t' -f 1,14 | tr '\t' '|' | cut -d'|' -f 1,6 | tr '|' '\t' > Milani_Sanger_GoodSequences.txt
```

It is important to double check the file to make sure that it is clean

```
AR049_4	k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pyronemataceae;g__Trichophaea;s__Trichophaea_sp
AR070_4	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales;f__Cantharellales_fam_Incertae_sedis;g__Sistotrema;s__Sistotrema_sp
AR065_4	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Inocybaceae;g__Inocybe;s__Inocybe_pseudorubens
AR076_4	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales;f__Clavulinaceae;g__Clavulina;s__Clavulina_sp
```

Now file have to be merged with the unite database

```
#Merging the fasta files
cat sh_refs_qiime_ver8_dynamic_s_04.02.2020.fasta ./Sanger_data/Milani_Sanger_GoodSequences.fasta > sh_refs_qiime_ver8_dynamic_s_04.02.2020_w_Sanger.fasta

#Merging taxonomy files
cat sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020.txt ./Sanger_data/Milani_Sanger_GoodSequences.txt > sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020_w_Sanger.txt
```

Now we can use this files for reassigning the OTUs including the sanger sequences. 

# Reassigning OTUs using novel database

## Assigning taxonomy

We can now align our OTU sequences to the UNITE database using Blast. 

```
parallel_assign_taxonomy_blast.py -O 16 -t /mnt/DataHDD/Data2/unite_qiime/sh_qiime_release_fungi_s_04.02.2020/sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020_w_Sanger.txt \
                                                                  -r /mnt/DataHDD/Data2/unite_qiime/sh_qiime_release_fungi_s_04.02.2020/sh_refs_qiime_ver8_dynamic_s_04.02.2020_w_Sanger.fasta \
                                                                  -i rep_seq.fasta \
                                                                  -o blast_re-assigned_taxonomy
```

Some OTUs were reassigned to sequences obtained from root tips
```
OTU_41	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Atheliales;f__Atheliaceae;g__Tylospora;s__Tylospora_sp_SH1648324.08FU	4e-122	AR115_4
OTU_27	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales;f__Cantharellales_fam_Incertae_sedis;g__Sistotrema;s__Sistotrema_sp_SH1539274.08FU	6e-149	AR015_4
OTU_48	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Thelephorales;f__Thelephoraceae;g__Pseudotomentella;s__Pseudotomentella_rhizopunctata_SH1638684.08FU	5e-17	AR078_4
OTU_51	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Atheliales;f__Atheliaceae;g__Tylospora;s__Tylospora_sp_SH1648324.08FU	3e-120	AR082_4
OTU_116	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Thelephorales;f__Thelephoraceae;g__Thelephora;s__Thelephora_terrestris_SH1502189.08FU	1e-165	AR012_4
OTU_98	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Boletales;f__Rhizopogonaceae;g__Rhizopogon;s__Rhizopogon_pseudoroseolus_SH1555177.08FU	4e-172	AR062_4
OTU_97	k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pyronemataceae;g__Trichophaea;s__Trichophaea_sp_SH1557658.08FU	2e-145	AR089_4
OTU_124	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Amanitaceae;g__Amanita;s__Amanita_muscaria_SH1553667.08FU	8e-167	AR118_4
OTU_142	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Atheliales;f__Atheliaceae;g__Tylospora;s__Tylospora_sp_SH1648324.08FU	3e-104	AR082_4
OTU_190	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Tricholomataceae;g__Tricholoma;s__Tricholoma_pessundatum_SH1681112.08FU	0.0	AR021_4
OTU_261	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales;f__Clavulinaceae;g__Clavulina;s__Clavulina_sp_SH1606319.08FU	2e-22	AR076_4
OTU_295	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Atheliales;f__Atheliaceae;g__Tylospora;s__Tylospora_sp_SH1648324.08FU	6e-112	AR082_4
OTU_292	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Tricholomataceae;g__Tricholoma;s__Tricholoma_pessundatum_SH1681112.08FU	4e-30	AR021_4
OTU_328	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Atheliales;f__Atheliaceae;g__Tylospora;s__Tylospora_sp_SH1648324.08FU	3e-107	AR082_4
OTU_392	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Atheliales;f__Atheliaceae;g__Tylospora;s__Tylospora_sp_SH1648324.08FU	6e-118	AR115_4
OTU_371	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Thelephorales;f__Thelephoraceae;g__Pseudotomentella;s__Pseudotomentella_rhizopunctata_SH1638684.08FU	6e-174	AR034_4
OTU_450	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales;f__Cantharellales_fam_Incertae_sedis;g__Sistotrema;s__Sistotrema_sp_SH1539274.08FU	1e-20	AR015_4
OTU_660	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Atheliales;f__Atheliaceae;g__Tylospora;s__Tylospora_sp_SH1648324.08FU	5e-134	AR115_4
OTU_764	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales;f__Cantharellales_fam_Incertae_sedis;g__Sistotrema;s__Sistotrema_sp_SH1539274.08FU	3e-19	AR015_4
OTU_899	k__Fungi;p__Ascomycota;c__Leotiomycetes;o__Helotiales;f__Myxotrichaceae;g__Oidiodendron;s__Oidiodendron_sp_SH1565480.08FU	4e-122	AR020_5
OTU_926	k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Amanitaceae;g__Amanita;s__Amanita_muscaria_SH1553667.08FU	2e-139	AR060_4
```

Based on these resutls, a new BIOM file will be made for analyzing the obtained OTUs from soil samples.

## Make otu table biom file

Using qiime1, a biom file is created and the taxonomy and sample data could be added:
````
#Create biom file
biom convert --table-type="OTU table" -i ../otu_table.txt -o otu_table_TomARv2.biom --to-json

#adding taxonomy
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy \
				  --observation-metadata-fp blast_re-assigned_taxonomy/rep_seq_tax_assignments.txt \
				  -i otu_table_TomARv2.biom -o otu_table_TomARv2.tax.biom --output-as-json
				  
#adding sample data
biom add-metadata -m read_map.txt -i otu_table_TomARv2.tax.biom -o otu_table_TomARv2.tax.sample.biom --output-as-json
```

## Filter NTC from samples
```
#create Biom for negative controls (NTCs)
filter_samples_from_otu_table.py -i otu_table_TomARv2.tax.sample.biom -o otu_table_NTC.biom -m read_map.txt -s "Samples:NTC"

#Biom without NTC samples
filter_samples_from_otu_table.py -i otu_table_TomARv2.tax.sample.biom -o otu_table_TomARv21.tax.sample.biom -m read_map.txt -s "Samples:AR.Tom"

#remove non zeros from NTCs
filter_otus_from_otu_table.py -i otu_table_NTC.biom -o otu_table_NTC_Non.biom -n 1

#Convert NTC_biom file to a tab-separated txt table
biom convert -i otu_table_NTC_Non.biom -o otus_to_remove.txt --to-tsv

#Filter NTC OTUS using txt file from NTC samples
filter_otus_from_otu_table.py -i otu_table_TomARv21.tax.sample.biom -o otu_table_TomasARv2-1.clean.biom -e otus_to_remove.txt

#Filter non-zero samples
filter_otus_from_otu_table.py -i otu_table_TomasARv2-1.clean.biom -o otu_table_TomasARv2-1.clean1.biom -n 1

#Check counts in final biom file
biom summarize-table -i otu_table_TomasARv2-1.clean1.biom
```

