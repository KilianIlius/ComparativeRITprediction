from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import SearchIO
import os


# Function to BLAST your query RITs against the genome they originated from:
def blast_xml(query, genome, output_xml_name, word_size, max_hsps):

    # compute appropriate e-value
    seq_count = 0
    default_e_value = 0.1

    for seq in SeqIO.parse(query, "fasta"):
        seq_count += 1

    e_value = default_e_value/seq_count

    # do the BLAST search
    cmd = "blastn -query " + query + " -db " + genome + " -outfmt 5 -out " + output_xml_name + " -word_size " + str(word_size) + " -evalue " + str(e_value) \
          + " -max_hsps " + str(max_hsps)

    os.system(cmd)


# Function to embed Query RITs within genomic sequence
def embed_rit(genome_name, rit_fa, extracted_rits_xml, out, emb_length):
    sbjct_start_list = []
    sbjct_end_list = []

    result_handle = open(extracted_rits_xml, "r")
    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sbjct_start_list.append(hsp.sbjct_start)
                sbjct_end_list.append(hsp.sbjct_end)
                break

    out_file = open(out, "w")

    i = 0

    for genome in SeqIO.parse(genome_name, "fasta"):
        for rit in SeqIO.parse(rit_fa, "fasta"):

            if sbjct_start_list[i] < sbjct_end_list[i]:

                result = genome.seq[sbjct_start_list[i] - emb_length:sbjct_end_list[i] + emb_length]
                out_file.write(str(">" + rit.description + "+\n" + result + "\n"))

                i += 1

            else:

                result = genome.seq[sbjct_end_list[i] - emb_length:sbjct_start_list[i] + emb_length]
                out_file.write(str(">" + rit.description + "-\n" + result + "\n"))

                i += 1

    out_file.close()


# Function to transmute a given sequence
def transmute(in_seq, out):
    out_file = open(out, "w")

    for seq in SeqIO.parse(in_seq, "fasta"):
        temp_seq = Seq(seq.seq)
        result_seq = ""

        for i in range(len(temp_seq)):
            if temp_seq[i] == "a" or temp_seq[i] == "A":
                result_seq += "G"
            elif temp_seq[i] == "c" or temp_seq[i] == "C":
                result_seq += "T"
            elif temp_seq[i] == "g" or temp_seq[i] == "G":
                result_seq += "C"
            elif temp_seq[i] == "t" or temp_seq[i] == "T":
                result_seq += "A"

        out_file.write(str(">" + seq.description + "\n" + result_seq + "\n"))

    out_file.close()


# Function to BLAST query against multiple genomes
def blast_rits_all_genomes(genomes, query_name, word_size, max_hsps, outfmt):
    # compute appropriate e-value
    seq_count = 0
    default_e_value = 0.1

    for seq in SeqIO.parse(query_name, "fasta"):
        seq_count += 1

    e_value = default_e_value / seq_count

    out_file_list = []

    for genome_name in genomes:
        out_filename = output_path + genome_name.replace("fasta", "") + "vs_rits.xml"

        cmd = "blastn -query " + query_name + " -db " + genome_name + " -out " + out_filename \
              + " -outfmt " + str(outfmt) + " -word_size " + str(word_size) + " -evalue " + str(e_value) + " -max_hsps " + str(max_hsps)

        os.system(cmd)
        out_file_list.append(out_filename)

    return out_file_list


# Function to extract and embed all hits found with BLAST
def extract_and_embed_blast_out(genome_list, blast_output_list, out_file_list, emb_length):
    k = 0

    for output_name in blast_output_list:
        sbjct_start_list = []
        sbjct_end_list = []
        query_name_list = []
        query_id = []
        i = 0
        j = 0
        genome = SeqIO.read(genome_list[k], "fasta")
        result_handle = open(output_name, "r")
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            i += 1
            for alignment in blast_record.alignments:

                for hsp in alignment.hsps:
                    sbjct_start_list.append(hsp.sbjct_start)
                    sbjct_end_list.append(hsp.sbjct_end)
                    query_id.append(i)

        for qresult in SearchIO.parse(output_name, "blast-xml"):
            j += 1
            for x in query_id:
                if x == j:
                    query_name_list.append(qresult.id)

        out_file = open(out_file_list[k], "w")

        for n in range(len(sbjct_start_list)):
            if sbjct_start_list[n] < sbjct_end_list[n]:
                out_file.write(str(">" + query_name_list[n] + "\n" + genome.seq[sbjct_start_list[n] - emb_length:
                                                                                sbjct_end_list[n] + emb_length] + "\n"))
            else:
                result = genome.seq[sbjct_end_list[n] - emb_length:sbjct_start_list[n] + emb_length]
                out_file.write(str(">" + query_name_list[n] + "\n" + result + "\n"))
        k += 1


# Define output path, path of input data, and filenames for output data
input_path = "test_files/"
output_path = "C:/Users/kilia/PycharmProjects/RNIE_score_averaging/test_files/example_transmutation_approach_out/"

input_rits = input_path + "test_rits_B._subtilis.fasta"
input_genome_name = "B._subtilis_NC_000964.3.fasta"
input_genome_path = input_path + input_genome_name
subject_genomes = ["B._subtilis_NC_000964.3.fasta", "B._tequilensis_NZ_CP048852.1.fasta", "B._vallismortis_NZ_CP033052.1.fasta"]
subject_genomes_trans = ["B._subtilis_NC_000964.3_trans.fasta", "B._tequilensis_NZ_CP048852.1_trans.fasta", "B._vallismortis_NZ_CP033052.1_trans.fasta"]

embedded_rits_out = output_path + "test_rits_B.sub_emb.fasta"
trans_rits_out = output_path + "test_rits_B.sub_emb_trans.fasta"

# Define length of sequence to be added downstream and upstream to the RIT/homolog sequence during embedding/BLAST extraction
embedding_length = 30


# BLAST search input RITs against Input genome and save the results in XML-format
rits_vs_genome_xml = output_path + "rits_vs_genome.xml"
blast_xml(input_rits, input_genome_name, rits_vs_genome_xml, 7, 1)

# Embed input RITs within genomic sequence
embed_rit(input_genome_path, input_rits, rits_vs_genome_xml, embedded_rits_out, embedding_length)

# Transmute embedded positives for negatives as well as all subject genomes
transmute(embedded_rits_out, trans_rits_out)

for item in subject_genomes:
    transmute(input_path + item, output_path + item.replace(".fasta", "") + "_trans.fasta")

# BLAST positives against normal genomes
blast_all_genomes_pos_result_list = blast_rits_all_genomes(subject_genomes, embedded_rits_out, 7, 1, 5)

# BLAST negatives against transmuted genomes
blast_all_genomes_neg_result_list = blast_rits_all_genomes(subject_genomes_trans, trans_rits_out, 7, 1, 5)

# Create lists of genome file paths
genomes_paths = []
extract_and_emb_pos_outfile_list = []
trans_genomes_paths = []
extract_and_emb_neg_outfile_list = []

for item in subject_genomes:
    genomes_paths.append(input_path + item)
    extract_and_emb_pos_outfile_list.append(output_path + item[:6] + "_rit_homologs.fasta")

for item in subject_genomes_trans:
    trans_genomes_paths.append(output_path + item)
    extract_and_emb_neg_outfile_list.append(output_path + item[:6] + "_neg_homologs.fasta")


# Extract and embed BLAST hits for positives
extract_and_embed_blast_out(genomes_paths, blast_all_genomes_pos_result_list, extract_and_emb_pos_outfile_list, embedding_length)

# Extract and embed BLAST hits for negatives
extract_and_embed_blast_out(trans_genomes_paths, blast_all_genomes_neg_result_list, extract_and_emb_neg_outfile_list, embedding_length)

