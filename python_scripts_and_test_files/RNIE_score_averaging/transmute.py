from Bio import SeqIO
from Bio.Seq import Seq


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


seq_to_transmute = "test_files/B._subtilis_NC_000964.3.fasta"
transmuted_seq_out = "test_files/B._subtilis_NC_000964.3_trans.fasta"

transmute(seq_to_transmute, transmuted_seq_out)
