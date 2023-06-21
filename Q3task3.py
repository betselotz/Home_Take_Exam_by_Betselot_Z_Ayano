import os
import Levenshtein

class LevenshteinDistanceCalculator:
    def __init__(self, file_path):
        self.file_path = file_path
        self.species_let7_sequences = {}

    def calculate_levenshtein_distance(self):
        if not os.path.isfile(self.file_path):
            print("The specified file does not exist.")
            return

        species_code = ""
        sequence = ""

        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    if species_code and sequence:
                        self.species_let7_sequences.setdefault(species_code, []).append(sequence)

                    species_code = line[1:].strip().split("let-7")[0]
                    sequence = ""
                else:
                    sequence = line.strip()

            if species_code and sequence:
                self.species_let7_sequences.setdefault(species_code, []).append(sequence)

        for species_code, sequences in self.species_let7_sequences.items():
            total_distance = 0
            total_sequences = len(sequences)

            if total_sequences < 2:
                continue

            for i in range(total_sequences - 1):
                for j in range(i + 1, total_sequences):
                    total_distance += Levenshtein.distance(sequences[i], sequences[j])

            average_distance = total_distance / (total_sequences * (total_sequences - 1) / 2)
            print(f"Average Levenshtein distance for let-7 miRNA in species '{species_code}': {average_distance:.2f}")

    def run_analysis(self):
        self.calculate_levenshtein_distance()


file_path = r"D:\1_MSC_AAU\2_Bioinformatics\1st_Year\2nd_Semester\4_Advanced_Programming_for_Bioinformatics\Betselot_Zerihun\Betselot_Zerihun\mature.fa"

calculator = LevenshteinDistanceCalculator(file_path)
calculator.run_analysis()
