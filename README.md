# EVEs annotation Refinement in reference genome

Pipeline developed by Elisa Pischedda, while in the Bonizzoni Lab at the University of Pavia (Italy).

## Purpose

The pipeline used in Whitfield et al. 2017, allow to find EVEs within reference genomes. Here, we extend the pipeline of Whitefield with a refinement of the annotation that allow the user to produce an accurate list of EVEs in a standard format.
Optionally the user can manually check the result and modify the scripts to obtaine the desired result.

### Reverse blastx
The input of the pipeline is based on a blastx obtained with one of the following commands. It is necessary to set the output file as described in the format parameter (-outfmt in blastx or -f in diamond):

```sh
blastx -query tophits.fasta \
-db NR/RefSeq_protein \
-evalue 1e-06 \
-outfmt '6 qseqid qstart qend salltitles saccver evalue qframe pident qcovs sstart send slen staxid' \
-out TopHits.blastx
```
```sh
diamond blastx -d nr_diamond.dmnd \
-e 1e-06 \
-f 6 qseqid qstart qend salltitles evalue qframe pident qcovhsp sstart send slen staxids \
-q TopHits.fasta -o TopHits.blastx &
```

- - - -

## Requirements
The pipeline uses the following programs.

* PYTHON 2.7 https://www.python.org/download/releases/
* PYTHON 3 https://www.python.org/download/releases/
* Virus-Host Classifier https://github.com/Kzra/VHost-Classifier
  It requires ETE3 toolkit http://etetoolkit.org/download/
  Note: the first time that it is used, it download the NCBI taxonomy database and put it in /home/username/.etetoolkit , thus the user have to install it in own account. To update the database use the following command:
```sh
python3 /ngs-data/Pipelines/Refine_EVE_annotation/Update_taxon_localdb.py

```

* check the latest version of Taxonkit https://bioinf.shenwei.me/taxonkit/
* update the NCBI taxonomy database for taxonkit with the following command and then uncompress. By default taxonkit requires the database to be in the "/home/username/.taxonkit" folder.

```sh
cd /home/username/.taxonkit
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
```

* download the file virushostdb.tsv from virus host database in VHost-Classifier-master directory

```sh
cd /ngs-data/Pipelines/Refine_EVE_annotation/VHost-Classifier-master
wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv
```

## Command line

```sh
bash Refine_EVE_Annotation.sh \
-pipeline_directory /pathTo/Refine_EVE_annotation_folder \
-reverse_blast_tool diamond \
--VHC_directory /pathTo/VHost-Classifier \
-file_blastx /pathTo/TopHits.blastx \
-file_bed_tophit /pathTo/TopHits.bed \
-output_directory /pathTo/Output_directory \
-taxonkit_exe /pathTo/Taxonkit_0.6/taxonkit
```

### Parameters

| parameter name                          | description                                                                             | default |
|-----------------------------------------|-----------------------------------------------------------------------------------------|---------|
| pipeline_directory                      | absolute path of the pipeline directory                                                 |         |
| tool                                    | reverse blastx tool used, options: blastx, diamond.                                          | blastx  |
| file_blastx                             | absolute path of the output file of the reverse blast                                   |         |
| file_bed_tophit                         | absolute path of the bed file of the tophits from the pipeline of Whitfield             |         |
| output_directory                        | absolute path of the output directory                                                   |         |
| taxonkit_exe                            | absolute path of blastn executable                                                       |         |

<br>

### Output tables fields

Two folders are created in the output path defined by the user. The ‘VHC’ folder contains the files created by VirusHost-Classifier while the ‘Output’ folder contains the tab separated value (tsv) tables produced by the pipeline. The ‘CompleteTable_classified’ file is the most informative and containes the following fields.

| Field                       | Description                                                                          |
|-----------------------------|--------------------------------------------------------------------------------------|
| ID                          | Unique ID of the EVE, assigned automatically                                         |
| Scaffold                    | Assembly scaffold/chromosome                                                         |
| Start                       | EVE start in the scaffold/chromosome                                                 |
| End                         | EVE end in the scaffold/chromosome                                                   |
| Len                         | EVE length (bp)                                                                      |
| q_seq_id                    | Original EVE ID (as in the user-provided BED file)                                   |
| Num_total_hit               | Total BLASTx hits                                                                    |
| Num_viral_hit               | Viral BLASTx hits                                                                    |
| Fraction_viral_hit          | Fraction of viral BLASTx hits on the total number of hits                            |
| Bestviral_Protein           | Best viral hit subject ID                                                            |
| Virus_name                  | Best viral hit subject viral species                                                 |
| BestviralEval               | Best viral hit E-Value                                                               |
| Num_subj_same_bestViralEval | Numbers of viral hits with the same E-Value                                          |
| Bestviral_frame             | Best viral hit frame                                                                 |
| Bestviral_percid            | Best viral hit percentage of identity (%)                                            |
| Bestviral_qcov              | Best viral hit query coverage (%)                                                    |
| Bestviral_start             | Best viral hit start on the subject                                                  |
| Bestviral_end               | Best viral hit end on the subject                                                    |
| Bestviral_slen              | Best viral hit total subject length                                                  |
| Bestviral_taxid             | Best viral hit subject NCBI taxonomy ID                                              |
| Bestviral_accession         | Best viral hit subject NCBI accession number                                         |
| BestSubj                    | Best non-viral hit subject ID                                                        |
| BestEval                    | Best non-viral hit subject E-Value                                                   |
| Num_subj_same_bestEval      | Numbers of non-viral hits with the same E-Value                                      |
| Sel_BestSubj                | Best non-viral hit from a non-predicted protein subject ID                           |
| Sel_BestEval                | Best non-viral hit from a non-predicted protein subject E-Value                      |
| Sel_Num_subj_same_bestEval  | Number of non-viral hit from a non-predicted protein subject with the   same E-Value |
| NIRVSmerged                 | Indicates if the EVE was merged from more overlapping hits                           |
| Viral Order                 | Best viral hit taxonomic order                                                       |
| Viral Family                | Best viral hit taxonmic family                                                       |
| Viral Order taxid           | NCBI taxonomy ID Best viral hit taxonomic order                                      |
| Viral Family taxid          | NCBI taxonomy ID Best viral hit taxonmic family                                      |

## References

[Whitfield et al., 2017]('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5698160' "Whitfield et al., 2017")

[Kitson et al., 2019]('https://academic.oup.com/bioinformatics/article-abstract/35/19/3867/5368011?redirectedFrom=fulltext' "Kitson et al., 2019")

<br>
