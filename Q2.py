class MiRNACounter:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7_count = 0

    def count_let7_miRNAs(self):
        with open(self.file_path, "r") as file:
            for line in file:
                if line.startswith(">"):
                    miRNA_name = line.strip()[1:]
                    if "let-7" in miRNA_name:
                        self.let7_count += 1

    def run_analysis(self):
        self.count_let7_miRNAs()
        print("Total number of let-7 miRNAs across all species:", self.let7_count)


file_path = r"D:\1_MSC_AAU\2_Bioinformatics\1st_Year\2nd_Semester\4_Advanced_Programming_for_Bioinformatics\Betselot_Zerihun\Betselot_Zerihun\mature.fa"

counter = MiRNACounter(file_path)
counter.run_analysis()
