def split_algn_file(algn_file, output_path):
    with open(algn_file, "r") as f:
        lines = f.read().split("//")

        f2 = open(algn_file, "r")
        all_lines = f2.readlines()

        sw_no_list = []

        for line in all_lines:
            if "ID" in line:
                sw_no_list.append(line.split("_")[2].strip())

        for i in range(len(lines) - 1):
            out_file = open(output_path + "/sw_" + str(sw_no_list[i]) + "_single_algn.sto", "w")
            out_file.write(lines[i] + "//")
            out_file.close()
            print(str(i) + "/" + str(len(lines) - 1))


multiple_algn_file = "//wsl.localhost/Ubuntu/home/kilian/hmmer_E._coli_sw/split1/s1_hmmer_algn.sto"
single_algn_output_path = "//wsl.localhost/Ubuntu/home/kilian/hmmer_E._coli_sw/split1/single_alignments/"

split_algn_file(multiple_algn_file, single_algn_output_path)
