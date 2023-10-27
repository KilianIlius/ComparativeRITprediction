from Bio.Seq import Seq
import os


def annotate_structure(rnie_infile, alignment_infile, out_path):

    # read RNIE cmsearch output file and store the results of seed-b in seed_b_lines
    with open(rnie_infile) as f:
        lines_cmsearch_file = f.readlines()
        seed_b_lines = []
        correct_seed = False

        for line in lines_cmsearch_file:
            if "seed-b" in line:
                correct_seed = True

            if correct_seed:
                seed_b_lines.append(line)

    # read alignment file save the window number in window_no and all lines in lines_align_file
    with open(alignment_infile) as g:
        window_no = 0
        lines_align_file = g.readlines()

        for line in lines_align_file:
            if "GF ID" in line:
                window_no = int(line.split()[2].split("_")[2])

    # find the best scoring cmsearch result for the correct sliding window in the seed-b results, save the orientation of the best hit and put out the corresponding
    # block of the cmsearch file including score, structure and aligned sequence
    k = 0
    scores_dict = {}

    for line in seed_b_lines:
        quit_window = False

        if quit_window:
            break

        if line == ">sliding_window_" + str(window_no) + "\n":
            k += 1
            while not quit_window:
                if "Plus" in seed_b_lines[k]:
                    k += 2
                    orientation = "+"
                    new_entry_exists = True

                    while new_entry_exists:
                        scores_dict[float(seed_b_lines[k + 1].split()[2].replace(",", ""))] = [orientation, k]

                        if "Query" not in seed_b_lines[k + 9]:
                            new_entry_exists = False
                            k -= 9

                        k += 9

                elif "Minus" in seed_b_lines[k]:
                    orientation = "-"
                    k += 2
                    new_entry_exists = True

                    while new_entry_exists:
                        scores_dict[float(seed_b_lines[k + 1].split()[2].replace(",", ""))] = [orientation, k]

                        if "Query" not in seed_b_lines[k + 9]:
                            new_entry_exists = False
                            quit_window = True
                            k -= 9

                        k += 9

                elif ">sliding_window_" in seed_b_lines[k]:
                    quit_window = True
                k += 1
        k += 1

    print(window_no)
    print(scores_dict)

    annotated_structure_alignment_outfile = open(out_path + "/sw_" + str(window_no) + "_single_algn_with_structure.sto", "w")

    genome_list = ("Escherichia", "Salmonella", "Shigella", "Citrobacter", "Kluyvera")

    if scores_dict == {}:
        print("no predicted RNIE structure found")
        first_line = False
        for line in lines_align_file:
            if line.startswith(genome_list) and not first_line:
                first_line = True
                algn_seq = line.strip().split()[1]

                white_space_line = ""
                line_parts = line.split(" ")

                for i in range(len(line_parts[0]) - 11):
                    white_space_line += " "

                for i in range(len(line_parts) - 2):
                    white_space_line += " "

                annotated_structure_alignment_outfile.write("#=GC SS_cons" + white_space_line)

                for i in range(len(algn_seq)):
                    annotated_structure_alignment_outfile.write(".")

                annotated_structure_alignment_outfile.write("\n")

            annotated_structure_alignment_outfile.write(line)
        return
    else:
        best_scoring_block = seed_b_lines[scores_dict[max(scores_dict)][1]:scores_dict[max(scores_dict)][1] + 7]

    # # add the structure from the cmsearch output to the alignment file

    # first get the pure structure and corresponding sequence without gaps
    predicted_structure = best_scoring_block[3].strip()
    predicted_structure_alignment = best_scoring_block[6].strip().split()[1]

    i = 0
    new_predicted_structure = ""
    new_predicted_structure_alignment = ""
    while i < len(predicted_structure_alignment):
        if predicted_structure_alignment[i] == "-":
            i += 1
        elif predicted_structure[i] == ":":
            i += 1
        else:
            new_predicted_structure += predicted_structure[i]
            new_predicted_structure_alignment += predicted_structure_alignment[i]
            i += 1

    new_predicted_structure_alignment = new_predicted_structure_alignment.upper()
    new_predicted_structure_alignment = Seq(new_predicted_structure_alignment.replace("U", "T"))
    new_predicted_structure = new_predicted_structure.replace(":", ".").replace("_", ".")

    if scores_dict[max(scores_dict)][0] == "-":
        new_predicted_structure_alignment = new_predicted_structure_alignment.reverse_complement()
        new_predicted_structure = new_predicted_structure[::-1].replace(">", "x").replace("<", ">").replace("x", "<")

    print(new_predicted_structure)
    print(new_predicted_structure_alignment)

    no_open = 0
    no_close = 0

    for char in new_predicted_structure:
        if char == "<":
            no_open += 1
        elif char == ">":
            no_close += 1

    temp_no = 0
    break_reached = False
    if no_open != no_close:
        if no_open > no_close:
            for char in new_predicted_structure:
                if char == "<" and not break_reached:
                    temp_no += 1
                    if temp_no == no_close:
                        break_reached = True
                elif char == "<" and break_reached:
                    new_predicted_structure = new_predicted_structure[:temp_no] + "." + new_predicted_structure[temp_no + 1:]
                    temp_no += 1
        elif no_open < no_close:
            for char in new_predicted_structure:
                if char == ">":
                    new_predicted_structure = new_predicted_structure[:temp_no] + "." + new_predicted_structure[temp_no + 1:]
                    no_close -= 1
                    if no_close == no_open:
                        break
                temp_no += 1

        print(new_predicted_structure)

    # add the structure

    first_line_detected = False
    for line in lines_align_file:
        if line.startswith(genome_list) and not first_line_detected:
            first_line_detected = True
            algn_seq = line.strip().split()[1]

            white_space_line = ""
            line_parts = line.split(" ")

            for i in range(len(line_parts[0]) - 11):
                white_space_line += " "

            for i in range(len(line_parts) - 2):
                white_space_line += " "

            annotated_structure_alignment_outfile.write("#=GC SS_cons" + white_space_line)

            algn_seq_pos = 0
            pred_struc_algn_pos = 0
            count = 0
            estm_end = (len(new_predicted_structure_alignment)/3)*2
            structure_found = False
            # print(algn_seq)

            while algn_seq_pos < len(algn_seq):
                while pred_struc_algn_pos < len(new_predicted_structure_alignment) and not structure_found:
                    # print(algn_seq_pos, pred_struc_algn_pos)
                    if algn_seq_pos >= len(algn_seq):
                        break
                    if algn_seq[algn_seq_pos] == ".":
                        algn_seq_pos += 1
                        if pred_struc_algn_pos != 0:
                            new_predicted_structure = new_predicted_structure[:pred_struc_algn_pos] + "x" + new_predicted_structure[pred_struc_algn_pos:]
                        else:
                            annotated_structure_alignment_outfile.write(".")
                    elif not algn_seq[algn_seq_pos] == new_predicted_structure_alignment[pred_struc_algn_pos]:
                        annotated_structure_alignment_outfile.write(".")
                        algn_seq_pos -= count - 1
                        pred_struc_algn_pos = 0
                        gap_no = new_predicted_structure.count("x")
                        new_predicted_structure = new_predicted_structure.replace("x", "")
                        annotated_structure_alignment_outfile.write(gap_no * ".")
                        count = 0
                    else:
                        pred_struc_algn_pos += 1
                        algn_seq_pos += 1
                        count += 1

                        if pred_struc_algn_pos >= estm_end:
                            if algn_seq_pos + len(new_predicted_structure) - estm_end < len(algn_seq):
                                annotated_structure_alignment_outfile.write(new_predicted_structure.replace("x", "."))
                            else:
                                annotated_structure_alignment_outfile.write(new_predicted_structure[0:(len(algn_seq) - algn_seq_pos + int(estm_end))].replace("x", "."))

                            structure_found = True
                            algn_seq_pos += (len(new_predicted_structure_alignment) - estm_end)

                if not algn_seq_pos == len(algn_seq):
                    annotated_structure_alignment_outfile.write(".")
                algn_seq_pos += 1

            annotated_structure_alignment_outfile.write("\n")

            annotated_structure_alignment_outfile.write(line)
        else:
            annotated_structure_alignment_outfile.write(line)


single_alignments_path = "splits/split1_gardner/nhmmer_alignments/single_alignments"
rnie_output_path = "splits/split1_gardner/split1_sliding_windows-geneMode.cmsearch"
alignments_with_structure_out_path = "splits/split1_gardner/nhmmer_alignments/single_alignments_with_structure"

for filename in os.listdir(single_alignments_path):
    fullname = os.path.join(single_alignments_path, filename)
    annotate_structure(rnie_output_path, fullname, alignments_with_structure_out_path)

