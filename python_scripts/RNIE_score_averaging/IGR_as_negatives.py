from Bio import SeqIO
from BCBio import GFF
from Bio.Blast import NCBIXML
from Bio import SearchIO
import os


# Function to extract all IGRs (with min length = min_len_igrs) not containing known terminators from a target genome and cut down IGRs containing terminators
# Output how many terminators were found in which IGR -> term_count_in_IGRs_calc.txt
def extract_igrs(genome, genome_annotation, rit_annotation, min_len_igrs, igrs_out, igrs_cut_out, out_path):

    # Create dictionary with all gene locations
    gene_loc_dict = {}
    sep = ":"
    k = 0

    in_file = genome_annotation
    in_handle = open(in_file)

    for rec in GFF.parse(in_handle):
        for feat in rec.features:
            clean_loc = str(feat.location)[1:-4]
            start_loc = int(clean_loc.split(sep, 1)[0])
            end_loc = int(clean_loc.split(sep, 1)[1])
            gene_loc_dict[k] = (start_loc, end_loc)
            k += 1

    # Create dictionary where overlapping genes are merged into one annotation
    gene_loc_dict_no_overlap = {}
    n = 0
    start_save = gene_loc_dict[0][0]
    end_save = gene_loc_dict[0][1]

    for i in range(len(gene_loc_dict) - 1):
        if gene_loc_dict[i][1] > gene_loc_dict[i + 1][0]:
            end_save = gene_loc_dict[i + 1][1]
        else:
            gene_loc_dict_no_overlap[n] = (start_save, end_save)
            start_save = gene_loc_dict[i + 1][0]
            end_save = gene_loc_dict[i + 1][1]
            n += 1

    # Extract IGRs from genome and delete IGRs containing termination sites
    content = []

    with open(rit_annotation) as f:
        for line in f:
            content.append(line.strip().split())

    for genome in SeqIO.parse(genome, "fasta"):
        pos = 0
        j = 0
        term_count = 0
        igr_len = 0
        min_len_igr = min_len_igrs
        igr_with_term_sites_free_space = []
        igr_lens = []
        out_file = open(igrs_out, "w")
        out_file_igr_count = open(out_path + "term_count_in_IGRs_calc.txt", "w")
        out_file_igr_positive_cut = open(igrs_cut_out, "w")

        for i in range(len(gene_loc_dict_no_overlap)):
            while j < len(content):
                if pos < int(content[j][1]) < gene_loc_dict_no_overlap[i][0] or pos < int(content[j][2]) < gene_loc_dict_no_overlap[i][0]:

                    if term_count != 0:
                        out_file_igr_count.write(str(term_count) + "\n")
                        igr_with_term_sites_free_space.append(igr_len - term_count * 100)

                    igr_len = gene_loc_dict_no_overlap[i][0] - pos
                    term_count = 1
                    out_file_igr_count.write("IGR_" + str(i) + " " + str(pos) + ":" + str(gene_loc_dict_no_overlap[i][0]) + "\n" +
                                             str(igr_len) + "\n")

                    # extract parts of IGRs containing termination sites (subtract 65 from each side of the termination site)
                    k = j
                    temp_pos = pos

                    while k < len(content):

                        if int(content[k][1]) >= gene_loc_dict_no_overlap[i][0]:
                            if temp_pos < gene_loc_dict_no_overlap[i][0] - min_len_igr:
                                out_file_igr_positive_cut.write("> cutIGR_" + str(temp_pos) + " " + str(gene_loc_dict_no_overlap[i][0]) + "\n"
                                                                + str(genome.seq[temp_pos:gene_loc_dict_no_overlap[i][0]]) + "\n")
                            break

                        if temp_pos < int(content[k][1]) - min_len_igr:
                            out_file_igr_positive_cut.write("> cutIGR_" + str(temp_pos) + " " + str(int(content[k][1])) + "\n"
                                                            + str(genome.seq[temp_pos:(int(content[k][1]))]) + "\n")
                            temp_pos = int(content[k][2])
                            k += 1
                        else:
                            temp_pos = int(content[k][2])
                            k += 1

                    pos = gene_loc_dict_no_overlap[i][1]
                    j += 1
                    break

                elif int(content[j][1]) > gene_loc_dict_no_overlap[i][0]:
                    if min_len_igr < len(genome.seq[pos:gene_loc_dict_no_overlap[i][0]]):
                        out_file.write("> IGR_" + str(i) + " " + str(pos) + ":" + str(gene_loc_dict_no_overlap[i][0]) +
                                       "\n" + str(genome.seq[pos:gene_loc_dict_no_overlap[i][0]]) + "\n")
                        igr_lens.append(len(genome.seq[pos:gene_loc_dict_no_overlap[i][0]]))
                    pos = gene_loc_dict_no_overlap[i][1]
                    break
                else:
                    j += 1
                    term_count += 1

        out_file.close()
        out_file_igr_positive_cut.close()
        out_file_igr_count.close()


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
def embed_rit(genome_name, rit_fa, extracted_rits_xml, out):
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

                result = genome.seq[sbjct_start_list[i] - 41:sbjct_end_list[i] + 40]
                out_file.write(str(">" + rit.description + "+\n" + result + "\n"))

                i += 1

            else:

                result = genome.seq[sbjct_end_list[i] - 41:sbjct_start_list[i] + 40]
                out_file.write(str(">" + rit.description + "-\n" + result + "\n"))

                i += 1

    out_file.close()


# Function to BLAST query against multiple genomes
def blast_rits_all_genomes(genomes, query_name, word_size, max_hsps, outfmt, split_name):
    # compute appropriate e-value
    seq_count = 0
    default_e_value = 0.1

    for seq in SeqIO.parse(query_name, "fasta"):
        seq_count += 1

    e_value = default_e_value / seq_count

    out_file_list = []

    for genome_name in genomes:
        out_filename = output_path + genome_name.replace(".fasta", "") + "_vs_" + split_name + ".xml"

        cmd = "blastn -query " + query_name + " -db " + genome_name + " -out " + out_filename \
              + " -outfmt " + str(outfmt) + " -word_size " + str(word_size) + " -evalue " + str(e_value) + " -max_hsps " + str(max_hsps)

        os.system(cmd)
        out_file_list.append(out_filename)

    return out_file_list


# Function to extract and embed all hits found with BLAST
def extract_and_embed_blast_out(genome_list, blast_output_list, out_file_list):
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
                if sbjct_start_list[n] - 30 < 1:
                    out_file.write(str(">" + query_name_list[n] + "\n" + genome.seq[1:sbjct_end_list[n] + 30] + "\n"))
                else:
                    out_file.write(str(">" + query_name_list[n] + "\n" + genome.seq[sbjct_start_list[n] - 30:
                                                                                    sbjct_end_list[n] + 30] + "\n"))
            else:
                result = genome.seq[sbjct_end_list[n] - 30:sbjct_start_list[n] + 30]
                out_file.write(str(">" + query_name_list[n] + "\n" + result + "\n"))
        k += 1


# Define output path, path of input data, and filenames for output data
input_path = "test_files/"
output_path = "test_files/example_IGR_as_negative_out/"

input_rits = input_path + "test_rits_B._subtilis.fasta"
input_genome_name = "B._subtilis_NC_000964.3.fasta"
input_genome_path = input_path + input_genome_name
input_genome_gene_annotation = "test_files/B._subtilis_genes.gff"
input_genome_rit_annotation = "test_files/RIT_annotations.bed"
subject_genomes = ["B._subtilis_NC_000964.3.fasta", "B._tequilensis_NZ_CP048852.1.fasta", "B._vallismortis_NZ_CP033052.1.fasta"]

# These negatives derived from IGRs serves as an example, in practice a split generated from all negatives should be used
igr_negative_split = "test_files/example_IGR_as_negative_out/example_negative_split.fasta"

# Output files of filtered IGRs not overlapping RITs
output_igrs_all = output_path + "B._sub_IGRs_no_RITs_min_len_50_complete.fasta"
embedded_rits_out = output_path + "test_rits_B.sub_emb.fasta"

# Define minimum length of negatives derived from filtered IGRs
igr_min_length = 50


# Extract IGRs not overlapping RITs with min length 50 to out_igrs; cut down IGRs containing RITs and save results with min length 50 to out_igrs_cut
out_igrs = output_path + "B._sub_IGRs_no_RITs_min_len_50.fasta"
out_igrs_cut = output_path + "B._sub_cutIGRs_min_len_50.fasta"
extract_igrs(input_genome_path, input_genome_gene_annotation, input_genome_rit_annotation, igr_min_length, out_igrs, out_igrs_cut, output_path)

# Merge files of IGRs_no_RITs_min_len_50 and cutIGRs_min_len_50
out_igrs_all = open(output_igrs_all, "w")
for line in open(out_igrs, "r"):
    out_igrs_all.write(line)

for line in open(out_igrs_cut, "r"):
    out_igrs_all.write(line)

# BLAST search input RITs against Input genome and save the results in XML-format
rits_vs_genome_xml = output_path + "rits_vs_genome.xml"
blast_xml(input_rits, input_genome_name, rits_vs_genome_xml, 7, 1)

# Embed input RITs within genomic sequence
embed_rit(input_genome_path, input_rits, rits_vs_genome_xml, embedded_rits_out)

# BLAST positives against all genomes
blast_all_genomes_pos_result_list = blast_rits_all_genomes(subject_genomes, embedded_rits_out, 7, 1, 5, "positives")

# BLAST negatives against all genomes
blast_all_genomes_neg_result_list = blast_rits_all_genomes(subject_genomes, igr_negative_split, 7, 1, 5, "negatives")


# Create lists of genome file paths
genomes_paths = []
extract_and_emb_pos_outfile_list = []
extract_and_emb_neg_outfile_list = []

for item in subject_genomes:
    genomes_paths.append(input_path + item)
    extract_and_emb_pos_outfile_list.append(output_path + item[:6] + "_rit_homologs.fasta")

for item in subject_genomes:
    extract_and_emb_neg_outfile_list.append(output_path + item[:6] + "_neg_homologs.fasta")


# Extract and embed BLAST hits for positives
extract_and_embed_blast_out(genomes_paths, blast_all_genomes_pos_result_list, extract_and_emb_pos_outfile_list)

# Extract and embed BLAST hits for negatives
extract_and_embed_blast_out(genomes_paths, blast_all_genomes_neg_result_list, extract_and_emb_neg_outfile_list)

