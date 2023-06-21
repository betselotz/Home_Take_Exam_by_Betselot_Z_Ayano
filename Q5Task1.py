import os
from Levenshtein import distance

class LevenshteinDistanceCalculator:
    def __init__(self, file_path):
        self.file_path = file_path
        self.let7_sequences = []
        self.total_let7_miRNA = 0
        self.total_distance = 0
        self.total_pairs = 0

    def calculate_levenshtein_distance(self):
        if not os.path.isfile(self.file_path):
            print("The specified file does not exist.")
            return

        with open(self.file_path, 'r') as file:
            species = ""
            sequence = ""
            for line in file:
                if line.startswith('>'):
                    header = line[1:].strip()
                    species = self.extract_species(header)
                    sequence = ""
                else:
                    sequence = line.strip()

                if 'let-7a' in header and species and sequence:
                    self.let7_sequences.append((species, sequence))
                    self.total_let7_miRNA += 1

        for i in range(len(self.let7_sequences) - 1):
            for j in range(i + 1, len(self.let7_sequences)):
                species1, seq1 = self.let7_sequences[i]
                species2, seq2 = self.let7_sequences[j]
                levenshtein_distance = distance(seq1, seq2)
                self.total_distance += levenshtein_distance
                self.total_pairs += 1
                print(f"Species: {species1} - {species2} | Levenshtein Distance: {levenshtein_distance}")

        if self.total_pairs > 0:
            average_distance = self.total_distance / self.total_pairs
            print(f"Total 'let-7a' miRNAs: {self.total_let7_miRNA}")
            print(f"Average Levenshtein Distance of 'let-7a' miRNAs: {average_distance:.2f}")

    def extract_species(self, header):
        species = header.split(' ')[0]
        return species



file_path = r"D:\1_MSC_AAU\2_Bioinformatics\1st_Year\2nd_Semester\4_Advanced_Programming_for_Bioinformatics\Betselot_Zerihun\Betselot_Zerihun\mature.fa"

calculator = LevenshteinDistanceCalculator(file_path)
calculator.calculate_levenshtein_distance()
