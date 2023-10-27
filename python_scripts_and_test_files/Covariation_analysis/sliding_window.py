from Bio import SeqIO


def sliding_window(sw_size, sw_step, igr_genome, real_genome_positions, windows_outfile):
    real_genome_pos = open(real_genome_positions, "r")
    pos_list = real_genome_pos.read().replace("[", "").replace("]", "").split(", ")
    real_genome_pos.close()

    out_file = open(windows_outfile, "w")

    for genome in SeqIO.parse(igr_genome, "fasta"):
        genome_len = len(genome.seq)
        k = 0
        start = 0
        end = sw_size
        step = sw_step

        while end < genome_len:
            out_file.write(">sliding_window_" + str(k) + " genome_start:" + str(pos_list[start]) + " genome_end:" + str(pos_list[end]) + " "
                           + str(pos_list[start:end]).replace("'", "").replace(" ", "") + "\n" + str(genome.seq[start:end]) + "\n")

            start += step
            end += step
            k += 1


window_size = 100
step_size = 20
merged_igr_split = "splits/split1_gardner/split1_merged_igrs.fasta"
IGRs_real_nucl_position_log = "splits/split1_gardner/split1_gardner_IGR_pos_translation.txt"
output_window_file = "splits/split1_gardner/split1_sliding_windows.fasta"


sliding_window(100, 20, merged_igr_split, IGRs_real_nucl_position_log, output_window_file)
