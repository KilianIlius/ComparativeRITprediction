# Prdicting Rho-independent terminators using comparative information

This GitHub page provides the data utilized and generated in the course of the associated Master Thesis. Two distinct approaches to classify RITs using comparative information were tested. Results of the RNIE score averaging approach can be viewed in the "Averaging_Approach" folder, while the results of the covariation analysis are listed in the "Covariation_Approach" folder. Furthermore, the relevant gene annotations and genomes downloaded from the NCBI website can be found in the folders "gene_annotation" and "genomes". The sets of RITs used as positives, provided by Gardner et al. (2011), S. Strobel (2020) and Chen et al. (2013), are localized in the folder "known_RITs". Splits of positive and negative data are listed in "splits".
Additionally, all python scripts developed during the study can be downloaded from "python_scripts_and_test_files". The following sections list the required software and python packages as well as a short documentation on how to employ the code on your data or on the data used in this thesis to replicate the results.

## Required software

To run RNIE, R-scape and HMMER a UNIX system is required. In this thesis the Windows-Subsystem for Linux Ubunutu was used, providing the 22.04.2 LTS release.

Software and version numbers:

* Python 3<br/>
* [BLAST-2.13.0+-x64-win64](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/)<br/>
<t/>BLAST was installed on a Windows system and integrated into python scripts via os module.<br/>
* [Infernal-1.0.2](http://eddylab.org/infernal/)<br/>
<t/>Installing Infernal version 1.0.2 is crucial to ensure RNIE to work properly. Using newer versions of Infernal might cause problems when utilizing RNIE. <br/>
* [Easel](http://eddylab.org/infernal/)<br/>
<t/>Integrated in Infernal-1.0.2.<br/>
* [HMMER-3.3.2](http://eddylab.org/software/hmmer/)<br/>
* [RNIE-0.01](https://github.com/ppgardne/RNIE)<br/>
* [R-scape-v2.0.0.q](http://eddylab.org/software/rscape/)<br/>
* [RNAfold](http://rna.tbi.univie.ac.at/)<br/>
<t/>Not necessary to run the scripts but useful to check the predicted MFE structure of RIT candidates

## Required Python 3 packages:<br/>

The following Python 3 packages are required to run the scripts:

* [itertools](https://docs.python.org/3/library/itertools.html)<br/>
* [matplotlib](https://matplotlib.org/)<br/>
* [numpy](https://numpy.org/)<br/>
* [os](https://docs.python.org/3/library/os.html)<br/>
* [biopython](https://biopython.org/)<br/>
* [bcbio-gff](https://pypi.org/project/bcbio-gff/)<br/>

## RNIE score averaging

The following serves as a quick tutorial on how to run the RNIE score averaging approach. This tutorial is employed on B. subtilis RITs. The tutorial files are described in brackets and can be found in the folder 'python_scripts_and_test_files\RNIE_score_averaging\test_files').<br/>
To start your RNIE score averaging approach you will need the following data:

* Split of known RITs as positives in FASTA format (test_rits_B._subtilis.fasta)<br/>
* Genome of the species your positives originated from in FASTA format (B._subtilis_NC_000964.3.fasta)<br/>
* Gene annotation for this genome in GFF format (B._subtilis_genes.gff)<br/>
* Genome annotation of the most complete set of RIT known for this species in BED format (RIT_annotations.bed)<br/>
* Genomes of species related to your start species in FASTA format (B._tequilensis_NZ_CP048852.1.fasta, B._vallismortis_NZ_CP033052.1.fasta)<br/>
<t/>In the study a set of 10 related species was choosen based on GC-content and RIT occurrence. For this tutorial two species related to B. subtilis suffice to ensure quick results.<br/>

To carry out later BLAST homology searches for each genome a BLAST database has to be created. (makeblastdb -in <genome.fasta> -dbtype nucl)


### Transmutation strategy

If you choose to use the transmutation strategy to generate negative data, you also need transmuted versions of all your genomes. Therefore, you can use <tt>transmute.py</tt>:

#### <tt>transmute.py</tt>

* seq_to_transmute: path to your input sequence<br/>
* transmuted_seq_out: define output filename<br/>

Generate BLAST databases for your transmuted genomes (makeblastdb -in <genome_trans.fasta> -dbtype nucl)

#### <tt>transmutation_approach.py</tt>

This script embeds your set of positives in genomic data, transmutes the results to receive a set of negatives and conducts a BLAST search to identify homologs for positives and negatives within all related species. Homologs identified by BLAST are then extracted and also embedded in genomic sequence.<br/>
The following files and parameters have to be set:

* input_path: path to your input data. (Genomes, positive set, etc.)<br/>
* output_path: path where the output data should be saved.<br/>
* input_rits: path to your positive RIT set.<br/>
* input_genome_name: name of the genome FASTA file your positives originated from.<br/>
* subject_genomes: list of genome FASTA files of your chosen species.<br/>
* subject_genomes_trans: list of output filenames for your transmuted genomes.<br/>
* embedded_rits_out: path and filename for the embedded positives output. <br/>
* trans_rits_out: path and filename for the generated negative set.<br/>
* embedding_length: Integer to define length of sequence to be added downstream and upstream of the RIT/homolog via embedding.<br/>

Homologs for positives and negatives are saved to your output directory. The homologs output files are named by the first 6 characters of the genome they were found in followed by '_rit_homologs' for positives and '_neg_homologs' for negatives. Also, the BLAST results for each genome are saved in XML format ('genome'_vs_positives/'genome'_vs_negatives). These BLAST outputs are mandatory for the later RNIE score averaging.

### IGRs as negatives approach

If you choose to use filtered IGRs as negatives you can use the <tt>IGR_as_negatives.py</tt> script to derive negatives and find/extract homologs for positives and negatives:

#### <tt>IGR_as_negatives.py</tt>

This script derives a set of negatives by filtering IGRs for sequences not overlapping with RIT annotations. A subset of these negatives is then used to carry out the homology search. In this example a small set of such IGR negatives is used (test_files/example_IGR_as_negative_out/example_negative_split.fasta). In practice a larger split of negatives should be used, generated by splitting the complete set of generated negatives (see 'splits\splits_B._sub_RITs\negatives_predicted_from_IGRs').<br/>
Most parameters to be set are equivalent to the transmutation approach, extra parameteres are listed below:

* input_genome_gene_annotation: path to your genome gene annotation file in GFF format.<br/>
* input_genome_rit_annotation: path to your genome RIT annotation file in BED format.<br/>
* igr_negative_split: split of negatives you want to use.<br/>
* output_igrs_all: path where the FASTA file containing all negatives derived from IGR filtering should be saved.<br/>
* igr_min_length: minimum length each negative should have.<br/>

Just as in the transmutation approach homologs found for positives and negatives within each genome are saved in FASTA format. Additionally the generated negatives are saved. Also, the BLAST results for each genome are saved in XML format ('genome'_vs_positives/'genome'_vs_negatives). These BLAST outputs are mandatory for the later RNIE score averaging.

### RNIE score computation and averaging

Regardless of whether you chose the transmutation approach or the IGR filtering approach, the XML files containing the BLAST results for positives and negatives against all genomes are needed for RNIE score averaging. RNIE score averaging can be employed using the <tt>RNIE_score_averaging.py</tt> script. The script must be run twice. In the first iteration the number of BLAST matches per sample is computed and plotted and two FASTA files are generated containing all homologs of positives/negatives. These FASTA files alongside the FASTA files of initial positive and negatives sequences are then given to RNIE to compute RNIE scores. The RNIE output files are then used in the second runthrough to compute the averaged RNIE scores, plot initial scores and averaged scores and create ROC curves depicting the calssification performance.

#### First iteration of <tt>RNIE_score_averaging.py</tt>

* workspace_path: directory for your output files<br/>
* path_positives: path to the BLAST XML files containing the homology search results for positives vs all genomes.<br/>
* path_negatives: path to the BLAST XML files containing the homology search results for negatives vs all genomes.<br/>
* path_genomes: path to a folder containing the genomes of all used species (also transmuted genomes if you chose the transmutation approach).<br/>
* out_file_positive_BLAST_matches: specify the filename for your output file containing all positve homologs.<br/>
* out_file_negative_BLAST_matches: specify the filename for your output file containing all negative homologs.<br/>

Use ./rnie.pl to compute RNIE scores for positive homologs, negative homologs, initial positive sequences, and initial negative sequences. Safe the GFF formated output of RNIE to your workspace.

#### Second iteration of <tt>RNIE_score_averaging.py</tt>

* ini_positives:<t/> path to your initial set of positive sequences.<br/>
* ini_negatives:<t/> path to your initial set of negative sequences.<br/>
* ini_scores_txt:<t/> serves as a log file to save all scores of initial sequences.<br/>
* avg_scores_txt:<t/> serves as a log file to save all averaged scores.<br/>
* rnie_out_pos:<t/> path to RNIEs GFF output of all positive homologs.<br/>
* rnie_out_neg:<t/> path to RNIEs GFF output of all negative homologs.<br/>
* rnie_out_initial_pos:<t/> path to RNIEs GFF output of all initial positive sequences.<br/>
* rnie_out_initial_neg:<t/> path to RNIEs GFF output of all initial negative sequences.<br/>

In the second iteration comment in lines 495-508 in <tt>RNIE_score_averaging.py</tt>. Run the script and you'll receive two histograms, depicting initial and averaged RNIE scores, and a ROC plot showing two ROC curves, one for classification using initial scores and one for classification using averaged scores. AUCs are plotted in the bottom right corner.

