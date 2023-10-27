from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from BCBio import GFF
import matplotlib.pyplot as plt
from Bio import SeqIO
from itertools import pairwise
import os


# Function to count BLAST matches found per positive/negative and plot; Create two fasta files containing all positive/negative homologs
def count_matches_and_create_rnie_file(path_pos, path_neg, genomes, out_pos, out_neg):
    blast_count_dict_positives = {}
    blast_count_dict_negatives = {}
    rnie_seqs_dict_positives = {}
    rnie_seqs_dict_negatives = {}

    # Iterate through xml files of positives and count BLAST matches; Create Dictionary with all hit sequences for RNIE
    for filename in os.listdir(path_pos):
        if not filename.endswith(".xml"):
            continue
        fullname = os.path.join(path_pos, filename)

        result_handle = open(fullname, "r")
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            query_name = blast_record.query
            query_name = query_name + "_"
            database = blast_record.database

            if query_name not in blast_count_dict_positives:
                blast_count_dict_positives[query_name] = 0

            if query_name not in rnie_seqs_dict_positives:
                rnie_seqs_dict_positives[query_name] = []

            k = 0
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    for genome in SeqIO.parse(genomes + database, "fasta"):
                        k += 1
                        if hsp.sbjct_start < hsp.sbjct_end:
                            if hsp.sbjct_start > 20:
                                rnie_seqs_dict_positives[query_name] = rnie_seqs_dict_positives[query_name] + \
                                                                       [str(genome.seq[int(hsp.sbjct_start) - 20:int(hsp.sbjct_end) + 20])]
                            else:
                                rnie_seqs_dict_positives[query_name] = rnie_seqs_dict_positives[query_name] + [str(genome.seq[0:int(hsp.sbjct_end) + 20])]
                        else:
                            if hsp.sbjct_end > 20:
                                sequence = Seq(genome.seq[int(hsp.sbjct_end) - 20:int(hsp.sbjct_start) + 20])
                            else:
                                sequence = Seq(genome.seq[0:int(hsp.sbjct_start) + 20])
                            sequence = sequence.reverse_complement()
                            rnie_seqs_dict_positives[query_name] = rnie_seqs_dict_positives[query_name] + [str(sequence)]
                blast_count_dict_positives[query_name] = blast_count_dict_positives[query_name] + k

    # Iterate through xml files of negatives and count BLAST matches; Create Dictionary with all hit sequences for RNIE
    for filename in os.listdir(path_neg):
        if not filename.endswith(".xml"):
            continue
        fullname = os.path.join(path_neg, filename)

        result_handle = open(fullname, "r")
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            query_name = blast_record.query
            query_name = query_name.strip().split()[0]
            database = blast_record.database

            if query_name not in blast_count_dict_negatives:
                blast_count_dict_negatives[query_name] = 0

            if query_name not in rnie_seqs_dict_negatives:
                rnie_seqs_dict_negatives[query_name] = []

            k = 0
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    for genome in SeqIO.parse(genomes + database, "fasta"):
                        k += 1
                        if hsp.sbjct_start < hsp.sbjct_end:
                            if hsp.sbjct_start > 20:
                                rnie_seqs_dict_negatives[query_name] = rnie_seqs_dict_negatives[query_name] + \
                                                                       [str(genome.seq[int(hsp.sbjct_start) - 20:int(hsp.sbjct_end) + 20])]
                            else:
                                rnie_seqs_dict_negatives[query_name] = rnie_seqs_dict_negatives[query_name] + [str(genome.seq[0:int(hsp.sbjct_end) + 20])]
                        else:
                            if hsp.sbjct_end > 20:
                                sequence = Seq(genome.seq[int(hsp.sbjct_end) - 20:int(hsp.sbjct_start) + 20])
                            else:
                                sequence = Seq(genome.seq[0:int(hsp.sbjct_start) + 20])
                            sequence = sequence.reverse_complement()
                            rnie_seqs_dict_negatives[query_name] = rnie_seqs_dict_negatives[query_name] + [str(sequence)]
                blast_count_dict_negatives[query_name] = blast_count_dict_negatives[query_name] + k

    # Compute average number of matches per sample and create a list containing the number of matches for each sample
    list_matches_positives = []
    list_matches_negatives = []
    sum_matches_positives = 0

    for key in blast_count_dict_positives:
        sum_matches_positives = sum_matches_positives + blast_count_dict_positives[key]
        list_matches_positives.append(blast_count_dict_positives[key])

    sum_matches_negatives = 0
    for key in blast_count_dict_negatives:
        sum_matches_negatives = sum_matches_negatives + blast_count_dict_negatives[key]
        list_matches_negatives.append(blast_count_dict_negatives[key])

    # Create a scatter-plot showing the number of matches for positives and negatives
    y_pos = range(len(list_matches_positives))
    y_neg = range(len(list_matches_negatives))
    plt.scatter(y_pos, list_matches_positives, c="b", marker=".", label="#matches per positive sample")
    plt.scatter(y_neg, list_matches_negatives, c="r", marker=".", label="#matches per negative sample")
    plt.legend(loc="upper right")
    plt.show()

    # Output all homologous sequences in a fasta file for RNIE (one for the positives, one for the negatives)

    out_pos_file = open(out_pos, "w")
    out_neg_file = open(out_neg, "w")

    for key in rnie_seqs_dict_positives:
        k = 0
        for value in rnie_seqs_dict_positives[key]:
            sequence = value.replace("-", "")
            out_pos_file.write("> " + str(key) + str(k) + "\n" + sequence + "\n")
            k += 1

    for key in rnie_seqs_dict_negatives:
        k = 0
        key_name = key.strip().split()
        for value in rnie_seqs_dict_negatives[key]:
            sequence = value.replace("-", "")
            out_neg_file.write("> " + str(key_name[0]) + "_" + str(k) + "\n" + sequence + "\n")
            k += 1

    out_neg_file.close()
    out_pos_file.close()


# Function to remove the notes from RNIEs GFF output files and returns file names
def remove_notes(file_name):
    sep = ";Note"

    out_file_name = file_name.replace(".gff", "") + "_nonotes.gff"
    out_file = open(out_file_name, "w")

    with open(file_name, 'r', encoding='utf-8') as file:
        data = file.readlines()

    for i in range(len(data)):
        data[i] = data[i].split(sep, 1)[0]
        out_file.write(data[i] + "\n")

    return out_file_name


# Function for RNIE score averaging
def rnie_score_averaging(rnie_pos, rnie_neg, rnie_pos_ini, rnie_neg_ini, pos_homologs, neg_homologs, ini_pos, ini_neg, out_ini, out_avg, bins, log=False):
    # Get scores from GFF file (RNIE-output)
    dict_rnie_scores_positives = {}
    dict_rnie_scores_negatives = {}

    in_file = rnie_pos
    in_handle = open(in_file)

    for rec in GFF.parse(in_handle):
        rec_id = rec.id
        best_score = 0

        for feat in rec.features:
            feat_score = str(feat.qualifiers["score"])[2:-2]
            feat_score = float(feat_score)

            if feat_score > best_score:
                best_score = feat_score

        dict_rnie_scores_positives[rec_id] = [best_score]

    in_file = rnie_neg
    in_handle = open(in_file)

    for rec in GFF.parse(in_handle):
        rec_id = rec.id
        best_score = -9

        for feat in rec.features:
            feat_score = str(feat.qualifiers["score"])[2:-2]
            feat_score = float(feat_score)

            if feat_score > best_score:
                best_score = feat_score

        dict_rnie_scores_negatives[rec_id] = [best_score]

    # fill dict_RNIE with 0s where no RNIE hit was found for a sequence
    sequence_ids_neg = []
    sequence_ids_pos = []

    for genome in SeqIO.parse(neg_homologs, "fasta"):
        sequence_ids_neg.append(genome.id)

    for genome in SeqIO.parse(pos_homologs, "fasta"):
        sequence_ids_pos.append(genome.id)

    for entry in sequence_ids_neg:
        if entry not in dict_rnie_scores_negatives.keys():
            dict_rnie_scores_negatives[entry] = [0]

    for entry in sequence_ids_pos:
        if entry not in dict_rnie_scores_positives.keys():
            dict_rnie_scores_positives[entry] = [0]

    all_scores_positives = []
    all_scores_negatives = []

    for key in dict_rnie_scores_positives:
        for value in dict_rnie_scores_positives[key]:
            all_scores_positives.append(value)

    for key in dict_rnie_scores_negatives:
        for value in dict_rnie_scores_negatives[key]:
            all_scores_negatives.append(value)

    # averaging RNIE scores (positives)
    average_scores_positives = []
    average_scores_positives_dict = {}

    for key in dict_rnie_scores_positives:
        key_name = str(key)
        key_name = key_name.split("_")[0]

        if key_name not in average_scores_positives_dict:
            average_scores_positives_dict[key_name] = []

        average_scores_positives_dict[key_name] = average_scores_positives_dict[key_name] + dict_rnie_scores_positives[key]

    for key in average_scores_positives_dict:
        k = 0
        temp_score = 0
        for value in average_scores_positives_dict[key]:
            temp_score = temp_score + value
            k += 1
        temp_score = temp_score/k
        average_scores_positives.append(temp_score)

    # averaging RNIE scores (negatives)
    average_scores_negatives = []
    average_scores_negatives_dict = {}

    for key in dict_rnie_scores_negatives:
        key_name = str(key)
        key_name = key_name.split("_")
        key_name = key_name[0] + "_" + key_name[1]

        if key_name not in average_scores_negatives_dict:
            average_scores_negatives_dict[key_name] = []

        average_scores_negatives_dict[key_name] = average_scores_negatives_dict[key_name] + dict_rnie_scores_negatives[key]

    for key in average_scores_negatives_dict:
        k = 0
        temp_score = 0
        for value in average_scores_negatives_dict[key]:
            temp_score = temp_score + value
            k += 1
        temp_score = temp_score/k
        average_scores_negatives.append(temp_score)

    # plot scores of initial sequences against averaged score

    dict_rnie_scores_ini_positives = {}
    dict_rnie_scores_ini_negatives = {}

    in_file = rnie_pos_ini
    in_handle = open(in_file)

    for rec in GFF.parse(in_handle):
        rec_id = rec.id
        best_score = 0

        for feat in rec.features:
            feat_score = str(feat.qualifiers["score"])[2:-2]
            feat_score = float(feat_score)

            if feat_score > best_score:
                best_score = feat_score

        dict_rnie_scores_ini_positives[rec_id] = best_score

    in_file = rnie_neg_ini
    in_handle = open(in_file)

    for rec in GFF.parse(in_handle):
        rec_id = rec.id
        best_score = 0

        for feat in rec.features:
            feat_score = str(feat.qualifiers["score"])[2:-2]
            feat_score = float(feat_score)

            if feat_score > best_score:
                best_score = feat_score

        dict_rnie_scores_ini_negatives[rec_id] = best_score

    pos_id_dict = {}
    neg_id_dict = {}

    for genome in SeqIO.parse(ini_pos, "fasta"):
        pos_id_dict[genome.id] = 0

    for genome in SeqIO.parse(ini_neg, "fasta"):
        neg_id_dict[genome.id] = 0

    for key in neg_id_dict:
        if key not in dict_rnie_scores_ini_negatives.keys():
            dict_rnie_scores_ini_negatives[key] = 0

    for key in pos_id_dict:
        if key not in dict_rnie_scores_ini_positives.keys():
            dict_rnie_scores_ini_positives[key] = 0

    ini_scores_pos_list = []
    ini_scores_neg_list = []

    # for key in average_scores_positives_dict:
    for key2 in dict_rnie_scores_ini_positives:
        ini_scores_pos_list.append(dict_rnie_scores_ini_positives[key2])

    # for key in average_scores_negatives_dict:
    for key2 in dict_rnie_scores_ini_negatives:
        ini_scores_neg_list.append(dict_rnie_scores_ini_negatives[key2])

    plt.hist(ini_scores_pos_list, bins=bins, label="initial positive scores", alpha=0.5, log=log)
    plt.hist(ini_scores_neg_list, bins=bins, label="initial negative scores", alpha=0.5, log=log)
    plt.legend(loc="upper right")
    plt.show()

    plt.hist(average_scores_positives, bins=bins, label="averaged RNIE scores positives", alpha=0.5, log=log)
    plt.hist(average_scores_negatives, bins=bins, label="averaged RNIE scores negatives", alpha=0.5, log=log)
    plt.legend(loc="upper right")
    plt.show()

    outfile_initial_scores = open(out_ini, "w")
    outfile_avg_scores = open(out_avg, "w")

    outfile_initial_scores.write("ini_pos_scores: \n")
    for item in ini_scores_pos_list:
        outfile_initial_scores.write(str(item) + "\n")
    outfile_initial_scores.write("ini_neg_scores: \n")
    for item in ini_scores_neg_list:
        outfile_initial_scores.write(str(item) + "\n")

    outfile_avg_scores.write("avg_pos_scores: \n")
    for item in average_scores_positives:
        outfile_avg_scores.write(str(item) + "\n")
    outfile_avg_scores.write("avg_neg_scores: \n")
    for item in average_scores_negatives:
        outfile_avg_scores.write(str(item) + "\n")

    outfile_initial_scores.close()
    outfile_avg_scores.close()

    return ini_scores_pos_list, ini_scores_neg_list, average_scores_positives, average_scores_negatives


# Function to compute two ROC curves (initial scores vs averaged scores)
def roc_curves(average_scores_positives, average_scores_negatives, ini_scores_pos_list, ini_scores_neg_list):

    def area_trapezoid(xs, ys):
        area = 0
        for (ax, ay), (bx, by) in pairwise(zip(xs, ys)):
            h = bx - ax
            area += h * by
        return area

    # Plot a ROC curve

    # choose which data to plot
    curve_scores_pos = average_scores_positives
    curve_scores_neg = average_scores_negatives

    num_positives = len(curve_scores_pos)
    num_negatives = len(curve_scores_neg)

    start_threshold = float(31)
    end_threshold = float(-100)
    step_size = float(0.1)
    score_range = start_threshold - end_threshold
    num_steps = int(score_range / step_size) + 2

    tpr_list = []
    fpr_list = []

    for i in range(num_steps):
        tp = 0
        fp = 0
        thresh = start_threshold - i * step_size

        for entry in curve_scores_pos:
            if entry > thresh:
                tp += 1
        tpr = tp / num_positives
        tpr_list.append(tpr)

        for entry in curve_scores_neg:
            if entry > thresh:
                fp += 1
        fpr = fp / num_negatives
        fpr_list.append(fpr)

    auc_avg = round(area_trapezoid(fpr_list, tpr_list), 4)
    legend_auc_avg = str("averaged scores; AUC = " + str(auc_avg))
    plt.plot(fpr_list, tpr_list, label=legend_auc_avg)

    # 2nd ROC curve in the same diagram
    curve_scores_pos = ini_scores_pos_list
    curve_scores_neg = ini_scores_neg_list

    num_positives = len(curve_scores_pos)
    num_negatives = len(curve_scores_neg)

    start_threshold = float(31)
    end_threshold = float(-100)
    step_size = float(0.1)
    score_range = start_threshold - end_threshold
    num_steps = int(score_range / step_size) + 2

    tpr_list = []
    fpr_list = []

    for i in range(num_steps):
        tp = 0
        fp = 0
        thresh = start_threshold - i * step_size

        for entry in curve_scores_pos:
            if entry > thresh:
                tp += 1
        tpr = tp / num_positives
        tpr_list.append(tpr)

        for entry in curve_scores_neg:
            if entry > thresh:
                fp += 1
        fpr = fp / num_negatives
        fpr_list.append(fpr)

    auc_ini = round(area_trapezoid(fpr_list, tpr_list), 4)
    legend_auc_ini = str("initial scores; AUC = " + str(auc_ini))
    plt.plot(fpr_list, tpr_list, label=legend_auc_ini)
    plt.legend()

    plt.show()


# specify path to workspace
workspace_path = "test_files/example_averaging/IGR_as_negatives_approach/"

# specify path to the BLAST results of positives/negatives against all genomes
path_positives = workspace_path + "BLAST_results/positives"
path_negatives = workspace_path + "BLAST_results/negatives"

# specify path to the genomes analyzed
path_genomes = "test_files/example_averaging/Genomes/"

# specify path to the initial sequences
ini_positives = workspace_path + "RNIE_results/initial_sequences/test_rits_B.sub_emb.fasta"
ini_negatives = workspace_path + "RNIE_results/initial_sequences/example_negative_split.fasta"

# specify output filenames
out_file_positive_BLAST_matches = workspace_path + "RNIE_results/homologs_positives.fasta"
out_file_negative_BLAST_matches = workspace_path + "RNIE_results/homologs_negatives.fasta"
ini_scores_txt = workspace_path + "initial_scores.txt"
avg_scores_txt = workspace_path + "average_scores.txt"

# Files after score computation by RNIE
rnie_out_pos = workspace_path + "RNIE_results/homologs_positives-geneMode-rnie.gff"
rnie_out_neg = workspace_path + "RNIE_results/homologs_negatives-geneMode-rnie.gff"
rnie_out_initial_pos = workspace_path + "RNIE_results/initial_sequences/test_rits_B.sub_emb-geneMode-rnie.gff"
rnie_out_initial_neg = workspace_path + "RNIE_results/initial_sequences/example_negative_split-geneMode-rnie.gff"


# Create two fasta files for RNIE containing all positive/negative homologous; Plot number of blast matches per sample
count_matches_and_create_rnie_file(path_positives, path_negatives, path_genomes, out_file_positive_BLAST_matches, out_file_negative_BLAST_matches)


# Use ./rnie.pl to compute RNIE scores for all positive and negative homologous.
# (DO NOT RUN THE SCRIPT BELOW BEFORE COMPUTING RNIE SCORES AND SPECIFYING FILE PATHS)


# # Remove notes from you RNIE output
# rnie_out_pos_no_notes = remove_notes(rnie_out_pos)
# rnie_out_neg_no_notes = remove_notes(rnie_out_neg)
# rnie_out_ini_pos_no_notes = remove_notes(rnie_out_initial_pos)
# rnie_out_ini_neg_no_notes = remove_notes(rnie_out_initial_neg)
#
# # Compute Average scores for positives and negatives; return list of initial scores and average scores for positives/negatives
# ini_scores_pos, ini_scores_neg, avg_scores_pos, avg_scores_neg = rnie_score_averaging(rnie_out_pos_no_notes, rnie_out_neg_no_notes,
#                                                                                       rnie_out_ini_pos_no_notes, rnie_out_ini_neg_no_notes,
#                                                                                       out_file_positive_BLAST_matches, out_file_negative_BLAST_matches,
#                                                                                       ini_positives, ini_negatives, ini_scores_txt, avg_scores_txt, 10)
#
# # Plot two ROC curves showing the classification with initial scores and with averaged scores
# roc_curves(avg_scores_pos, avg_scores_neg, ini_scores_pos, ini_scores_neg)
