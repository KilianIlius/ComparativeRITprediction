from Bio import SeqIO
from BCBio import GFF


# Function to extract all IGRs of a given genome, IGR file is named after the input Genome. Gene annotation file must be provided
def extract_igr(genome_file, gene_file, out_dir):
    # determine genome length
    for genome in SeqIO.parse(genome_file, "fasta"):
        genome_len = len(genome.seq)

    # Create dictionary with all gene locations
    gene_loc_dict = {}
    sep = ":"
    k = 0

    in_file = gene_file
    in_handle = open(in_file)

    for rec in GFF.parse(in_handle):
        for feat in rec.features:
            clean_loc = str(feat.location)[1:-4]
            start_loc = int(clean_loc.split(sep, 1)[0])
            end_loc = int(clean_loc.split(sep, 1)[1])
            gene_loc_dict[k] = (start_loc, end_loc)
            k += 1
        gene_loc_dict[k] = [genome_len, genome_len + 1]

    # Create dictionary where overlapping genes are merged into one annotation
    gene_loc_dict_no_overlap = {}
    n = 0
    start_save = gene_loc_dict[0][0]
    end_save = gene_loc_dict[0][1]

    for i in range(len(gene_loc_dict)):
        if i+1 not in gene_loc_dict.keys():
            gene_loc_dict[i+1] = [genome_len+2, genome_len+3]

        if gene_loc_dict[i][1] >= gene_loc_dict[i+1][0]:
            end_save = gene_loc_dict[i+1][1]
        else:
            gene_loc_dict_no_overlap[n] = (start_save, end_save)
            start_save = gene_loc_dict[i+1][0]
            end_save = gene_loc_dict[i+1][1]
            n += 1

    # Extract IGRs from genome

    for genome in SeqIO.parse(genome_file, "fasta"):
        pos = 0
        j = 0
        igr_lens = []
        genus = genome_file.strip().split("/")[1].split("_")[0] + "_" + genome_file.strip().split("/")[1].split("_")[1]
        out_file = out_dir + genus + "_IGRs.fasta"
        out_file_igrs_fasta = open(out_file, "w")

        for i in range(len(gene_loc_dict_no_overlap)):
            if pos < gene_loc_dict_no_overlap[i][0]:

                igr_lens.append(gene_loc_dict_no_overlap[i][0] - pos)

                out_file_igrs_fasta.write(">IGR_" + genus + "_" + str(j) + " start:" + str(pos) + " end:" + str(gene_loc_dict_no_overlap[i][0] - 1) + "\n" +
                                          str(genome.seq[pos:gene_loc_dict_no_overlap[i][0]]) + "\n")

                pos = gene_loc_dict_no_overlap[i][1]

                j += 1

        out_file_igrs_fasta.close()

    return out_file


# Function to concatenate given IGRs
def concatenate_igrs(igr_fasta, merged_igr_output, out_path):
    out_file_merged_igr = open(merged_igr_output, "w")

    genome_name = igr_fasta.strip().split("/")[1]
    genome_name = genome_name.strip().split("_")[0] + "_" + genome_name.strip().split("_")[1]

    out_file_real_genome_nt_pos = open(out_path + genome_name + "_IGR_pos_translation.txt", "w")

    merged_seq = ""
    real_genome_nt_pos_list = []

    out_file_merged_igr.write(">" + genome_name + "_IGR_genome" + "\n")

    for genome in SeqIO.parse(igr_fasta, "fasta"):

        merged_seq = merged_seq + genome.seq

        nt_start = genome.description.strip().split()[1]
        nt_start = int(nt_start.replace("start:", ""))

        nt_end = genome.description.strip().split()[2]
        nt_end = int(nt_end.replace("end:", ""))

        for num in range(nt_start, nt_end + 1):
            real_genome_nt_pos_list.append(num)

    out_file_merged_igr.write(str(merged_seq))
    out_file_real_genome_nt_pos.write(str(real_genome_nt_pos_list))
    out_file_real_genome_nt_pos.close()
    out_file_merged_igr.close()


# Input files; Important: genome_fasta must be named in the following way: "Genus_species_...fasta"
genome_fasta = "Entero_genomes/Escherichia_coli_genome.fasta"
gene_annotation = "Entero_gene_annotation/Escherichia_coli_gene_annotation.gff3"

igrs_output_path = "Entero_IGRs/"

# Output file
merged_igrs = "Merged_IGRs/Escherichia_coli_merged_IGR.fasta"
merged_igrs_out_path = "Merged_IGRs/"


igrs_out = extract_igr(genome_fasta, gene_annotation, igrs_output_path)

concatenate_igrs(igrs_out, merged_igrs, merged_igrs_out_path)
