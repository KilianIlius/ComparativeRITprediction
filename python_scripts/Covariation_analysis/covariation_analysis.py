from Bio import SeqIO
from Bio.Blast import NCBIXML
import os
import matplotlib.pyplot as plt
from itertools import pairwise


# Function to assign each window to either be positive, negative, or to be discarded
def assign_windows(positives_vs_genome_xml, windows_fasta, all_rits_vs_all_windows_xml, out_dir):
    positive_windows_list = []
    negative_windows_list = []
    discard_windows_list = []
    positives_positions_dict = {}

    result_handle = open(positives_vs_genome_xml, "r")
    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
        query_name = blast_record.query
        hsp_list = []

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.sbjct_start < hsp.sbjct_end:
                    temp_tuple = (hsp.sbjct_start, hsp.sbjct_end)
                else:
                    temp_tuple = (hsp.sbjct_end, hsp.sbjct_start)
                hsp_list.append(temp_tuple)

        positives_positions_dict[query_name] = hsp_list

    for window in SeqIO.parse(windows_fasta, "fasta"):
        window_name = window.description.strip().split()[0]
        window_start = int(window.description.strip().split()[1].replace("genome_start:", ""))
        window_end = int(window.description.strip().split()[2].replace("genome_end:", ""))
        window_nt = window.description.strip().split()[3].replace("[", "").replace("]", "").split(",")
        window_nt = [eval(i) for i in window_nt]

        for pos_key in positives_positions_dict:
            for pos_tuple in positives_positions_dict[pos_key]:
                if window_start < pos_tuple[0] < window_end < pos_tuple[1] or pos_tuple[0] < window_start < pos_tuple[1] < window_end:
                    if window_name not in positive_windows_list and window_name not in discard_windows_list:
                        discard_windows_list.append(window_name)
                    break
                elif pos_tuple[0] in window_nt and pos_tuple[1] in window_nt:
                    if window_name not in discard_windows_list and window_name not in positive_windows_list:
                        positive_windows_list.append(window_name)
                    break

        if window_name not in positive_windows_list and window_name not in discard_windows_list:
            negative_windows_list.append(window_name)

    j = 0
    new_positive_windows_list = []
    temp_list = []

    while j < len(positive_windows_list):
        k = j
        n = 0

        while k < len(positive_windows_list):
            if int(positive_windows_list[k].split("_")[2]) == int(positive_windows_list[j].split("_")[2]) + n:
                temp_list.append(positive_windows_list[k])
                k += 1
                n += 1
            else:
                # print(temp_list)
                middle_entry = (len(temp_list)) // 2
                new_positive_windows_list.append(temp_list[middle_entry])
                temp_list = []
                break

        j += n

    # for entry in negative_windows_list:
    #     print(entry)

    # delete windows from negative list where BLAST finds a hit for a terminator

    list_windows_containing_rits = []
    new_negative_windows_list = []

    result_handle2 = open(all_rits_vs_all_windows_xml, "r")
    blast_records2 = NCBIXML.parse(result_handle2)

    for blast_record2 in blast_records2:
        for alignment2 in blast_record2.alignments:
            list_windows_containing_rits.append(str(alignment2.hit_def).strip().split()[0])

    for entry in negative_windows_list:
        if entry not in list_windows_containing_rits:
            new_negative_windows_list.append(entry)

    # print lists of positive and negative windows
    out_pos_windows = open(out_dir + "pos_windows.txt", "w")
    out_neg_windows = open(out_dir + "neg_windows.txt", "w")

    for pos_window in new_positive_windows_list:
        out_pos_windows.write(str(pos_window) + "\n")

    for neg_window in new_negative_windows_list:
        out_neg_windows.write(str(neg_window) + "\n")

    out_pos_windows.close()
    out_neg_windows.close()

    return new_positive_windows_list, new_negative_windows_list


# Assign covariation e-values to each window
def get_cov_score(pos_window_list, neg_window_list, r_scape_path):
    pos_scores_dict = {}
    neg_scores_dict = {}
    no_pair_score = float(1.01)

    for filename in os.listdir(r_scape_path):
        if filename.endswith("sorted.cov"):
            fullname = os.path.join(r_scape_path, filename)
            window_num_r_scape = filename.strip().split("_")[2].split(".")[0]

            for window in pos_window_list:
                window_num_pos_list = window.strip().split("_")[2]

                if window_num_r_scape == window_num_pos_list:

                    with open(fullname) as f:
                        lines = f.readlines()

                    for line in lines:
                        if not line.startswith("#"):
                            if not line.strip() == "no significant pairs":
                                pos_scores_dict[window_num_pos_list] = float(line.strip().split()[4])
                                break
                            else:
                                pos_scores_dict[window_num_pos_list] = no_pair_score

            for window in neg_window_list:
                window_num_neg_list = window.strip().split("_")[2]

                if window_num_r_scape == window_num_neg_list:

                    with open(fullname) as f:
                        lines = f.readlines()

                    for line in lines:
                        if not line.startswith("#"):
                            if not line.strip() == "no significant pairs":
                                neg_scores_dict[window_num_neg_list] = float(line.strip().split()[4])
                                # print(window)
                                break
                            else:
                                neg_scores_dict[window_num_neg_list] = no_pair_score

    # print(pos_scores_dict)
    # print(neg_scores_dict)

    pos_scores_list = []
    neg_scores_list = []

    pos_no_cov_ratio = 0
    neg_no_cov_ratio = 0

    for key in pos_scores_dict:
        if not pos_scores_dict[key] == no_pair_score:
            pos_scores_list.append(pos_scores_dict[key])
        else:
            pos_scores_list.append(pos_scores_dict[key])
            pos_no_cov_ratio += 1

    pos_no_cov_ratio = pos_no_cov_ratio / len(pos_scores_dict)

    for key in neg_scores_dict:
        if not neg_scores_dict[key] == no_pair_score:
            neg_scores_list.append(neg_scores_dict[key])
        else:
            neg_scores_list.append(neg_scores_dict[key])
            neg_no_cov_ratio += 1

    neg_no_cov_ratio = neg_no_cov_ratio / len(neg_scores_dict)

    avg_pos_scores = sum(pos_scores_list)/len(pos_scores_list)
    avg_neg_scores = sum(neg_scores_list)/len(neg_scores_list)

    plt.hist(pos_scores_list, bins=50, label="positives \nno covariance ratio: " + str("%.3f" % pos_no_cov_ratio) + "\n" + "average e-value: " +
             str("%.3f" % avg_pos_scores), alpha=0.5, log=True)
    plt.hist(neg_scores_list, bins=50, label="negatives \nno covariance ratio: " + str("%.3f" % neg_no_cov_ratio) + "\n" + "average e-value: " +
             str("%.3f" % avg_neg_scores), alpha=0.5, log=True)
    plt.legend(loc="upper right")
    plt.show()

    sig_count_pos = 0
    for item in pos_scores_list:
        if item <= 0.05:
            sig_count_pos += 1

    print("Proportion of positives with e-values less than 0.05: " + str(sig_count_pos / len(pos_scores_list)))

    sig_count_neg = 0
    for item in neg_scores_list:
        if item <= 0.05:
            sig_count_neg += 1

    print("Proportion of negatives with e-values less than 0.05: " + str(sig_count_neg / len(neg_scores_list)))

    return pos_scores_list, neg_scores_list


# calculate AUC
def area_trapezoid(xs, ys):
    area = 0
    for (ax, ay), (bx, by) in pairwise(zip(xs, ys)):
        h = bx - ax
        area += h * by
    return "%.2f" % area


# Plot a ROC curve
def roc_curve(curve_scores_pos, curve_scores_neg):
    num_positives = len(curve_scores_pos)
    num_negatives = len(curve_scores_neg)

    start_threshold = float(1.02)
    end_threshold = float(0)
    step_size = float(0.01)
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

    plt.plot(fpr_list, tpr_list, label="e-values of nhmmer alignments; AUC: " + str(area_trapezoid(fpr_list, tpr_list)))
    plt.legend()
    plt.show()


# Define input files
split_path = "splits/split1_gardner/"

split_rits_vs_genome_xml = split_path + "BLAST_results/split1_vs_E_coli.xml"
all_rits_vs_genome_xml = split_path + "BLAST_results/all_e_coli_rits_vs_split1_windows.xml"

split_windows = split_path + "split1_sliding_windows.fasta"

r_scape_output_path = "splits/split1_gardner/r-scape_output"


# Determine positive and negative windows
split1_pos_window_list, split1_neg_window_list = assign_windows(split_rits_vs_genome_xml, split_windows, all_rits_vs_genome_xml, split_path)

# Assign covariation e-values to each window and plot histograms
pos_scores1, neg_scores1 = get_cov_score(split1_pos_window_list, split1_neg_window_list, r_scape_output_path)

# Plot ROC curve
roc_curve(pos_scores1, neg_scores1)
