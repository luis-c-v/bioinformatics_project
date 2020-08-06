"""
Name: Luis Cano Vazquez
Course: CSCI 1101B
Assignment: Project 4 - Bioinformatics
Date: 05/5/2020

Description: Takes Bioinformation from file, converts to a DNA to RNA strand and returns protein.
Bugs: None identified.
"""


CODON_DICT = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
              'UGC': 'C', 'UGU': 'C',
              'GAC': 'D', 'GAU': 'D',
              'GAA': 'E', 'GAG': 'E',
              'UUC': 'F', 'UUU': 'F',
              'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
              'CAC': 'H', 'CAU': 'H',
              'AUA': 'I', 'AUC': 'I', 'AUU': 'I',
              'AAA': 'K', 'AAG': 'K',
              'UUA': 'L', 'UUG': 'L', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
              'AUG': 'M',
              'AAC': 'N', 'AAU': 'N',
              'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
              'CAA': 'Q', 'CAG': 'Q',
              'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
              'AGC': 'S', 'AGU': 'S', 'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
              'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
              'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
              'UGG': 'W',
              'UAC': 'Y', 'UAU': 'Y'}


class DNA:
    """Creates an object that only analyzes DNA"""
    
    def __init__(self, strand):
        """creates a new DNA object"""
        self.strand = strand.lower()
        
    def validate(self):
        """checks whether the letters given to the DNA object are valid nucleotides"""
        nucl = 'acgt'
        for letter in self.strand:
            if letter not in nucl:
                only_nucl = False
                break
            else:
                only_nucl = True
        return only_nucl
    
    def __str__(self):
        """returns a human legiable line of infomation"""
        return self.strand.upper()
    
    def convert_to_rna(self):
        """returns a string but with the Ts replaced by Us"""
        return self.strand.replace('t', 'u').upper()
        
    def reverse_complement(self):
        """return a string that contains the reverse complement of the DNA sequence"""
        i = -1
        new_series = ''
        while abs(i) <= len(self.strand):
            # adds a reversed complement to the str new_series instead of attempting
            # to replace letters and indexes (which would cause confusion)
            if self.strand[i] == 'a':
                new_series += 't'
            elif self.strand[i] == 't':
                new_series += 'a'
            elif self.strand[i] == 'g':
                new_series += 'c'
            elif self.strand[i] == 'c':
                new_series += 'g'
            i -= 1
        return new_series.upper()


class RNA:
    """Creates an object that only takes in RNA"""
    
    def __init__(self, strand):
        """creates a new RNA object"""
        self.strand = strand.lower()
        
    def validate(self):
        """checks whether the letters given to the DNA object are valid nucleotides"""
        nucl = 'acgu'
        for letter in self.strand:
            if letter not in nucl:
                only_nucl = False
                # once the loop does find a non-nucleotide it changes T/F statement and breaks loop
                break
            else:
                only_nucl = True
        return only_nucl
    
    def __str__(self):
        """returns a human legiable line of infomation"""
        return self.strand.upper()
    
    def to_protein(self, index):
        """translates an RNA sequence into a protein"""
        self.strand = self.strand[index:].upper()
        # splits the string to specified index so we always start at '0'
        protein = ''     
        i = 0
        while i <= len(self.strand):
            tlc = self.strand[i:i + 3]
            # pulls only a codon at a time
            stop = ['UAG', 'UGA', 'UAA']
            if tlc not in stop:
                protein += CODON_DICT.get(tlc, '')
                #ensures that TLC is within start and stop codon
            else:
                break
            i += 3
        return protein


def read_dna_file(filename):
    """process each line of DNA data and print out the resulting protein str"""
    try:
        infile = open(filename)
        list_lines = infile.readlines()
        for dna in list_lines:
            dna = dna.strip().split()
            if dna != []:
                # pulls valid DNA strings
                test = DNA(dna[1])
                if test.validate():
                    # ensures that str is a DNA strand
                    rna = test.convert_to_rna()
                    new_rna = RNA(rna)
                    print(new_rna.to_protein(int(dna[0])))
                    # converts rna to protein using second item in list of DNA information
                else:
                    print('Invalid sequence: ' + dna[1])
    except FileNotFoundError:
        print('File not found.')


read_dna_file("DNASequences.txt")
read_dna_file("not_a_real_file.txt")