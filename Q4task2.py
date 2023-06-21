import Levenshtein
import matplotlib.pyplot as plt

class LevenshteinDistanceAnalyzer:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7a_mirnas = []
        self.distances = []

    def extract_let7a_mirnas(self):
        with open(self.file_path, "r") as file:
            miRNA_name = ""
            sequence = ""
            for line in file:
                if line.startswith(">"):
                    if miRNA_name != "":
                        if miRNA_name.startswith("hsa-let-7"):
                            self.let7a_mirnas.append((miRNA_name, sequence))
                    miRNA_name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()

            if miRNA_name.startswith("-let-7"):
                self.let7a_mirnas.append((miRNA_name, sequence))

    def calculate_pairwise_distances(self):
        for mirna, seq in self.let7a_mirnas:
            print("Pairwise Levenshtein distances for", mirna)
            for other_mirna, other_seq in self.let7a_mirnas:
                if mirna != other_mirna:
                    distance = Levenshtein.distance(seq, other_seq)
                    self.distances.append(distance)
                    print(f"Levenshtein distance between {mirna} and {other_mirna}: {distance}")

    def generate_histogram_plot(self):
        plt.hist(self.distances, bins=10)  # Adjust the number of bins as needed
        plt.xlabel("Levenshtein Distance")
        plt.ylabel("Frequency")
        plt.title("Distribution of Levenshtein Distances")
        plt.show()

    def run_analysis(self):
        self.extract_let7a_mirnas()
        self.calculate_pairwise_distances()
        self.generate_histogram_plot()


file_path = r"D:\1_MSC_AAU\2_Bioinformatics\1st_Year\2nd_Semester\4_Advanced_Programming_for_Bioinformatics\Betselot_Zerihun\Betselot_Zerihun\mature.fa"

analyzer = LevenshteinDistanceAnalyzer(file_path)
analyzer.run_analysis()
