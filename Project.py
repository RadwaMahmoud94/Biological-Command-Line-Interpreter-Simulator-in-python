import sys
import getopt
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


def cal_nbases(seq):
    c = 0
    for i in seq:
        if i == 'N' or i == 'n':
            c = c + 1
    return c
def is_valid(seq, type):
    if type == "DNA":
        dna = ['A', 'G', 'C', 'T']
        for i in seq:
            if i.upper() not in dna:
                return False
        return True
    elif type == "PROTEIN":
        protein = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        for i in seq:
            if i.upper() not in protein:
                return False
        return True
    elif type == "RNA":
        rna = ['A', 'G', 'C', 'U']
        for i in seq:
            if i.upper() not in rna:
                return False
        return True
    else:
        return "The type you entered is not right !"
def filter_nbases(seq):
    seq = seq.tomutable()
    for i in seq:
        if i == 'N' or i == 'n':
            seq.remove(i)
    return seq
def merge_fasta(files, output):
    with open(output, "w") as f:
        for file in files:
            with open(file) as infile:
                f.write(infile.read())
            f.write("\n")
        f.close()
    return f
def convert_to_fasta(genbank):
    with open(genbank) as input_handle, open("genbank_file.fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        SeqIO.write(sequences, output_handle, "fasta")
        input_handle.close()
        output_handle.close()
def online_alignment(seq):
    file = ""
    result_handle = NCBIWWW.qblast("blastn", "nt", seq)
    blast_record = NCBIXML.read(result_handle)
    e_value_thresh = 0.01
    for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                    if hsp.expect < e_value_thresh:
                        file += '****Alignment****\n'
                        file += 'sequence:' + str(alignment.title) + '\n'
                        file += 'length:' + str(alignment.length) + '\n'
                        file += 'e value:' + str(hsp.expect) + '\n'
                        file += str(hsp.query) + '\n'
                        file += str(hsp.match) + '\n'
                        file += str(hsp.sbjct) + '\n'
            file += '\n'
    return file


argv = sys.argv[1:]
opts, args = getopt.gnu_getopt(argv, 'o:')

commands = ["gc", "transcribe", "reverse_complement", "cal_nbases", "is_valid", "filter_nbases", "seq_alignment",
            "seq_alignment_files", "merge_fasta", "convert_to_fasta", "online_alignment"]
cmd = args[0]
args.remove(cmd)

if cmd in commands:
    if len(opts) == 0:
        opts = [('a', 'b')]
    for opt, arg in opts:
        try:
            if cmd == "gc":
                print("GC = ", GC(Seq(args[0])))


            elif cmd == "transcribe":
                print("Transcribe = ", (Seq(args[0])).transcribe())


            elif cmd == "reverse_complement":
                print("Reverse Complement =  ", Seq(args[0]).reverse_complement())


            elif cmd == "cal_nbases":
                print("Number of N-Bases = ", cal_nbases(Seq(args[0])))


            elif cmd == "filter_nbases":
                print("After Filteration = ", filter_nbases(Seq(args[0])))

            elif cmd == "is_valid":
                print(is_valid(Seq(args[0]), args[1].upper()))


            elif cmd == "seq_alignment":
                seq1 = Seq(args[0])
                seq2 = Seq(args[1])
                result = pairwise2.align.globalxx(seq1, seq2)
                if opt in ['-o']:
                    with open(arg, 'w') as f:
                        f.write(str(result))
                        f.close()
                else:
                    print(format_alignment(*result[0]))


            elif cmd == "seq_alignment_files":
                seq1 = SeqIO.read(args[0], "fasta")
                seq2 = SeqIO.read(args[1], "fasta")
                result = pairwise2.align.globalxx(seq1, seq2)
                if opt in ['-o']:
                    with open(arg, 'w') as f:
                        f.write(str(result))
                        f.close()
                else:
                    print(format_alignment(*result[0]))


            elif cmd == "merge_fasta":
                if opt in ['-o']:
                    merge_fasta(args, arg)
                else:
                    for file in args:
                        seq = SeqIO.parse(file, 'fasta')
                        for s in seq:
                            print(s.id)
                            print(s.seq)


            elif cmd == "online_alignment":
                seq = Seq(args[0])
                align = online_alignment(seq)

                if opt in ['-o']:
                    with open(arg, 'w') as f:
                        f.write(align)
                        f.close()
                else:
                    print(align)


            elif cmd == "convert_to_fasta":
                file_genbank = args[0]
                convert_to_fasta(file_genbank)

        except:
            print("Something went wrong !\n"
                  "please check command parameter(s)")
else:
    print("The command you entered is not found !")