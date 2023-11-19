import tkinter as tk
from tkinter import messagebox
import re

dna_codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

rna_codon_table = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

def valid_sequence(sequence, seq_type):
    valid_nucleotides = "ATGC" if seq_type == "1" else "AUGC"
    return all(nucleotide in valid_nucleotides for nucleotide in sequence)

def find_start_codon(frame, sequence):
    start_codons = ["ATG", "AUG"]
    for i in range(frame, len(sequence), 3):
        codon = sequence[i:i + 3]
        if codon in start_codons:
            return i
    return None

def translate(sequence, seq_type, codon_table):
    sequence = sequence.upper()
    translated_frames = []

    for frame in range(3):
        start_index = find_start_codon(frame, sequence)
        if start_index is not None:
            translated_frames.append([])
            for i in range(start_index, len(sequence), 3):
                codon = sequence[i:i + 3]
                if len(codon) == 3 and codon in codon_table:
                    if codon_table[codon] == "*":
                        break
                    translated_frames[frame].append(codon_table[codon])
                else:
                    translated_frames[frame].append("-")
        else:
            translated_frames.append(["No start codon found" if all(codon not in sequence for codon in ["ATG", "AUG"]) else "Start codon not in frame"])

    reverse_complement = sequence[::-1].replace("U", "T").translate(str.maketrans("ATGC", "TACG"))

    for frame in range(3):
        start_index = find_start_codon(frame, reverse_complement)
        if start_index is not None:
            translated_frames.append([])
            for i in range(start_index, len(reverse_complement), 3):
                codon = reverse_complement[i:i + 3]
                if len(codon) == 3 and codon in codon_table:
                    if codon_table[codon] == "*":
                        break
                    translated_frames[frame + 3].append(codon_table[codon])
                else:
                    translated_frames[frame + 3].append("-")
        else:
            translated_frames.append(["No start codon found" if all(codon not in reverse_complement for codon in ["ATG", "AUG"]) else "Start codon not in frame"])

    return translated_frames

class TranslationGUI:
    def __init__(self):
        self.window = tk.Tk()
        self.window.title("DNA/RNA Translator")
        self.window.configure(bg="#FFCFFF")
        self.window.geometry("730x630")
        self.window.resizable(False, False)
        self.seq_type_var = tk.StringVar(value="1")

        self.heading = tk.Label(self.window, text="DNA/RNA Translator", font=("Arial", 18, "bold"), fg="purple", pady=20, bg="#FFCFFF")
        self.heading.grid(row=0, column=0, columnspan=3)

        self.heading_label = tk.Label(self.window, text="Select Sequence Type:", font=("Arial", 14, "bold"), fg="purple", pady=20, padx=20, bg="#FFCFFF")
        self.heading_label.grid(row=1, column=0, columnspan=1)

        self.dna_button = tk.Radiobutton(self.window, text="DNA", variable=self.seq_type_var, value="1", bg="#FFCFFF")
        self.rna_button = tk.Radiobutton(self.window, text="RNA", variable=self.seq_type_var, value="2", bg="#FFCFFF")

        self.dna_button.grid(row=1, column=1)
        self.rna_button.grid(row=1, column=2)

        self.sequence_label = tk.Label(self.window, text="Sequence:", font=("Arial", 14, "bold"), fg="purple", pady=20, bg="#FFCFFF")
        self.sequence_entry = tk.Text(self.window, height=10, width=50, padx=5, pady=5)

        self.translate_button = tk.Button(self.window, text="Translate", command=self.translate, font=("Arial", 12, "bold"), bg="green", fg="white")

        self.output_label = tk.Label(self.window, text="ORF Translations:", font=("Arial", 14, "bold"), fg="purple", bg="#FFCFFF")
        self.output_text = tk.Text(self.window, height=10, width=50, padx=5, pady=5)

        self.clear_all_button = tk.Button(self.window, text="Clear Fields", command=self.clear_fields, font=("Arial", 12, "bold"), bg="red", fg="white")

        self.layout()

    def layout(self):
        self.sequence_label.grid(row=2, column=0)
        self.sequence_entry.grid(row=2, column=1, columnspan=2)
        
        sequence_vscrollbar = tk.Scrollbar(self.window, command=self.sequence_entry.yview)
        sequence_vscrollbar.grid(row=2, column=3, sticky='ns')
        self.sequence_entry.config(yscrollcommand=sequence_vscrollbar.set)

        self.translate_button.grid(row=3, column=1, pady=20)

        self.output_label.grid(row=4, column=0)
        self.output_text.grid(row=4, column=1, columnspan=2)

        output_vscrollbar = tk.Scrollbar(self.window, command=self.output_text.yview)
        output_vscrollbar.grid(row=4, column=3, sticky='ns')
        self.output_text.config(yscrollcommand=output_vscrollbar.set)

        self.clear_all_button.grid(row=5, column=1, pady=20)

    def clean_sequence(self, sequence):
        seq_type = self.seq_type_var.get()
        valid_characters = "ATGC" if seq_type == "1" else "AUGC"
        return re.sub(f"[^{valid_characters}]+", "", sequence)

    def translate(self):
        seq_type = self.seq_type_var.get()
        sequence = self.sequence_entry.get("1.0", tk.END).strip()
        sequence = self.clean_sequence(sequence)

        if seq_type == "1" and "ATGC" not in sequence:
            messagebox.showerror("Invalid Sequence", "The sequence you entered is not a valid DNA sequence! Please enter correct sequence.")
            return
        elif seq_type == "2" and "AUGC" not in sequence:
            messagebox.showerror("Invalid Sequence", "The sequence you entered is not a valid RNA sequence! Please enter correct sequence.")
            return

        codon_table = dna_codon_table if seq_type == "1" else rna_codon_table
        
        translation_result = translate(sequence, seq_type, codon_table)

        self.output_text.delete("1.0", tk.END)
        for i, frame in enumerate(translation_result, 1):
            direction = "5' to 3'" if i <= 3 else "3' to 5'"
            frame_number = i if i <= 3 else i - 3
            self.output_text.insert(
                tk.END, f"Frame {frame_number} ({direction}): {''.join(frame)}\n\n"
            )

    def clear_fields(self):
        self.sequence_entry.delete("1.0", tk.END)
        self.output_text.delete("1.0", tk.END)

if __name__ == "__main__":
    gui = TranslationGUI()
    gui.window.mainloop()
