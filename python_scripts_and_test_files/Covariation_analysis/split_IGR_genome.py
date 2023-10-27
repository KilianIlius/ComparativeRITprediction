from Bio import SeqIO
from Bio.Blast import NCBIXML
from generate_IGR_genome import concatenate_igrs


# Function to assign positive RITs to the IGRs they occur in and split them into two not overlapping groups
def find_positive_igrs(pos_vs_igrs_xml):
    igr_list = []
    split1_dict = {}
    split2_dict = {}
    out_file = open("log.txt", "w")
    contains = False
    m = 1
    del_list = []

    result_handle = open(pos_vs_igrs_xml, "r")
    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
        query_name = blast_record.query
        for alignment in blast_record.alignments:
            igr_list.append(alignment.hit_def)

        if m % 2 == 1:

            for item in igr_list:
                for key in split2_dict.keys():
                    for value in split2_dict[key]:
                        if item == value:
                            contains = True
                            break

            if not contains:
                split1_dict[query_name] = igr_list
                m += 1
            else:
                split2_dict[query_name] = igr_list
                for item in igr_list:
                    for key in split1_dict.keys():
                        for value in split1_dict[key]:
                            if item == value:
                                split2_dict[key] = split1_dict[key]
                                del_list.append(key)

                for entry in del_list:
                    if entry in split1_dict.keys():
                        del split1_dict[entry]

                del_list = []

        else:

            for item in igr_list:
                for key in split1_dict.keys():
                    for value in split1_dict[key]:
                        if item == value:
                            contains = True
                            break

            if not contains:
                split2_dict[query_name] = igr_list
                m += 1
            else:
                split1_dict[query_name] = igr_list
                for item in igr_list:
                    for key in split2_dict.keys():
                        for value in split2_dict[key]:
                            if item == value:
                                split1_dict[key] = split2_dict[key]
                                del_list.append(key)

                for entry in del_list:
                    if entry in split2_dict.keys():
                        del split2_dict[entry]

                del_list = []

        igr_list = []
        contains = False

    for key in split1_dict:
        for value in split1_dict[key]:
            for key2 in split2_dict:
                for value2 in split2_dict[key2]:
                    if value == value2:
                        print("Same value occurs in both splits")
                        print(value)

    out_file.write(str(split1_dict) + "\n\n\n" + str(split2_dict))

    return split1_dict, split2_dict


# Function to extract two splits of positives and their IGRs from the dictionaries computed with function find_positive_igrs
def extract_positive_splits(split1_dict, split2_dict, igr_file, rit_file, out_file_split1_rits, out_file_split1_igrs, out_file_split2_rits, out_file_split2_igrs):
    out_s1_rits = open(out_file_split1_rits, "w")
    out_s1_igrs = open(out_file_split1_igrs, "w")
    out_s2_rits = open(out_file_split2_rits, "w")
    out_s2_igrs = open(out_file_split2_igrs, "w")

    check_doubles_s1 = []
    check_doubles_s2 = []

    for rit in SeqIO.parse(rit_file, "fasta"):
        if rit.description in split1_dict.keys():
            out_s1_rits.write(">" + str(rit.description) + "\n" + str(rit.seq) + "\n")

        elif rit.description in split2_dict.keys():
            out_s2_rits.write(">" + str(rit.description) + "\n" + str(rit.seq) + "\n")

    for igr in SeqIO.parse(igr_file, "fasta"):
        for key in split1_dict:
            for value in split1_dict[key]:
                if igr.description == value and igr.description not in check_doubles_s1:
                    out_s1_igrs.write(">" + str(igr.description) + "\n" + str(igr.seq) + "\n")
                    check_doubles_s1.append(igr.description)

        for key2 in split2_dict:
            for value2 in split2_dict[key2]:
                if igr.description == value2 and igr.description not in check_doubles_s2:
                    out_s2_igrs.write(">" + str(igr.description) + "\n" + str(igr.seq) + "\n")
                    check_doubles_s2.append(igr.description)

    all_igrs = check_doubles_s1 + check_doubles_s2

    out_s1_rits.close()
    out_s1_igrs.close()
    out_s2_rits.close()
    out_s2_igrs.close()

    return all_igrs


# Function to extract all negative IGRs
def extract_negatives(igr_ids, igr_file, out_file_neg_igrs):
    out_neg_igr = open(out_file_neg_igrs, "w")

    for igr in SeqIO.parse(igr_file, "fasta"):
        if igr.description not in igr_ids:
            out_neg_igr.write(">" + str(igr.description) + "\n" + str(igr.seq) + "\n")


# Function to split negative IGRs into four splits
def split_negatives(infile_all_negatives, split1_pos, split2_pos, split3_pos, split4_pos, out_file_split1_neg, out_file_split2_neg,
                    out_file_split3_neg, out_file_split4_neg):
    split1_len = 0
    split2_len = 0
    split3_len = 0
    split4_len = 0

    j = 1

    out1 = open(out_file_split1_neg, "w")
    out2 = open(out_file_split2_neg, "w")
    out3 = open(out_file_split3_neg, "w")
    out4 = open(out_file_split4_neg, "w")

    for igr in SeqIO.parse(split1_pos, "fasta"):
        split1_len += len(igr.seq)

    for igr in SeqIO.parse(split2_pos, "fasta"):
        split2_len += len(igr.seq)

    for igr in SeqIO.parse(split3_pos, "fasta"):
        split3_len += len(igr.seq)

    for igr in SeqIO.parse(split4_pos, "fasta"):
        split4_len += len(igr.seq)

    for igr in SeqIO.parse(infile_all_negatives, "fasta"):
        if j % 4 == 1:
            out1.write(">" + str(igr.description) + "\n" + str(igr.seq) + "\n")
            split1_len += len(igr.seq)
            j += 1

        elif j % 4 == 2:
            out2.write(">" + str(igr.description) + "\n" + str(igr.seq) + "\n")
            split2_len += len(igr.seq)
            j += 1

        elif j % 4 == 3:
            out3.write(">" + str(igr.description) + "\n" + str(igr.seq) + "\n")
            split3_len += len(igr.seq)
            j += 1

        elif j % 4 == 0:
            out4.write(">" + str(igr.description) + "\n" + str(igr.seq) + "\n")
            split4_len += len(igr.seq)
            j += 1

    out1.close()
    out2.close()
    out3.close()
    out4.close()


# Function to pair positives and negatives into four splits and sort each split by IGR occurrence in the genome
def merge_pos_and_neg_igrs(pos_igrs, neg_igrs, out):
    out_file = open(out, "w")
    all_igrs_dict = {}

    for pos_igr in SeqIO.parse(pos_igrs, "fasta"):
        all_igrs_dict[pos_igr.description] = pos_igr.seq

    for neg_igr in SeqIO.parse(neg_igrs, "fasta"):
        all_igrs_dict[neg_igr.description] = neg_igr.seq

    dict_keys = list(all_igrs_dict.keys())
    dict_keys.sort(key=lambda x: int(x.split(" ")[0].split("_")[3]))

    sorted_all_igrs_dict = {i: all_igrs_dict[i] for i in dict_keys}

    for key in sorted_all_igrs_dict:
        out_file.write(">" + str(key) + "\n" + str(sorted_all_igrs_dict[key]) + "\n")

    out_file.close()


# INPUT
# Define paths to the known RITs
rits_gardner = "E._coli_rits/col_rits_gardner.fasta"
rits_chen = "E._coli_rits/col_rits_chen.fasta"
# Define path to BLAST XML results of each set of known RITs against the target 'IGR genome'
blast_gardner_xml = "BLAST_results/RITs_vs_E._coli_IGRs/e._coli_gardner_RITs_vs_e._coli_IGRs.xml"
blast_chen_xml = "BLAST_results/RITs_vs_E._coli_IGRs/e._coli_chen_RITs_vs_e._coli_IGRs.xml"

# Define path to the target genome IGRs generated with generate_IGR_genome.py
e_coli_igrs = "Entero_IGRs/Escherichia_coli_IGRs.fasta"

# OUTPUT
# Define paths to each split
split1_path = "splits/split1_gardner/"
split2_path = "splits/split2_gardner/"
split3_path = "splits/split3_chen/"
split4_path = "splits/split4_chen/"
# Define paths for each split of positive RITs
split1_rits = split1_path + "split1_rits_g.fasta"
split2_rits = split2_path + "split2_rits_g.fasta"
split3_rits = split3_path + "split3_rits_c.fasta"
split4_rits = split4_path + "split4_rits_c.fasta"
# Define paths for each split of positive IGRs
split1_pos_igrs = split1_path + "split1_igrs_g_positives.fasta"
split2_pos_igrs = split2_path + "split2_igrs_g_positives.fasta"
split3_pos_igrs = split3_path + "split3_igrs_c_positives.fasta"
split4_pos_igrs = split4_path + "split4_igrs_c_positives.fasta"
# Define paths for each split of negative IGRs
split1_neg_igrs = split1_path + "split1_igrs_negatives.fasta"
split2_neg_igrs = split2_path + "split2_igrs_negatives.fasta"
split3_neg_igrs = split3_path + "split3_igrs_negatives.fasta"
split4_neg_igrs = split4_path + "split4_igrs_negatives.fasta"
# Define paths for each split of positive/negative IGRs
split1_all_igrs = split1_path + "split1_igrs.fasta"
split2_all_igrs = split2_path + "split2_igrs.fasta"
split3_all_igrs = split3_path + "split3_igrs.fasta"
split4_all_igrs = split4_path + "split4_igrs.fasta"
# Define paths for each split of concatenated IGRs
split1_conc_igrs = split1_path + "split1_merged_igrs.fasta"
split2_conc_igrs = split2_path + "split2_merged_igrs.fasta"
split3_conc_igrs = split3_path + "split3_merged_igrs.fasta"
split4_conc_igrs = split4_path + "split4_merged_igrs.fasta"
# Define path for all negative igrs
all_neg = "splits/all_negative_igrs.fasta"


# Identify IGRs containing positives and split them into independent groups
split1_gardner_dict, split2_gardner_dict = find_positive_igrs(blast_gardner_xml)
split1_chen_dict, split2_chen_dict = find_positive_igrs(blast_chen_xml)

# Extract positive splits and save them in FASTA format (save IGRs and RIT sequences)
igrs_g = extract_positive_splits(split1_gardner_dict, split2_gardner_dict, e_coli_igrs, rits_gardner, split1_rits, split1_pos_igrs, split2_rits, split2_pos_igrs)
igrs_c = extract_positive_splits(split1_chen_dict, split2_chen_dict, e_coli_igrs, rits_chen, split3_rits, split3_pos_igrs, split4_rits, split4_pos_igrs)

# Save all positive IGR IDs
all_igr_ids = igrs_g + igrs_c

# Extract all IGRs not containing RITs as negatives
extract_negatives(all_igr_ids, e_coli_igrs, all_neg)

# Split negatives into four independent groups
split_negatives(all_neg, split1_pos_igrs, split2_pos_igrs, split3_pos_igrs, split4_pos_igrs, split1_neg_igrs, split2_neg_igrs, split3_neg_igrs, split4_neg_igrs)

# Pair each positive split with its negative split
merge_pos_and_neg_igrs(split1_pos_igrs, split1_neg_igrs, split1_all_igrs)
merge_pos_and_neg_igrs(split2_pos_igrs, split2_neg_igrs, split2_all_igrs)
merge_pos_and_neg_igrs(split3_pos_igrs, split3_neg_igrs, split3_all_igrs)
merge_pos_and_neg_igrs(split4_pos_igrs, split4_neg_igrs, split4_all_igrs)

# Concatenate each splits IGRs
concatenate_igrs(split1_all_igrs, split1_conc_igrs, split1_path)
concatenate_igrs(split2_all_igrs, split2_conc_igrs, split2_path)
concatenate_igrs(split3_all_igrs, split3_conc_igrs, split3_path)
concatenate_igrs(split4_all_igrs, split4_conc_igrs, split4_path)
