# bioinformatics_project


# Bioinformatics involves the use of (and development of) computer science tools and techniques to acquire, store, 
# organize, archive, analyze and visualize biological data. 


# **** TASK 1 ****
# Create a DNA class to model a strand of DNA.
# Your class should be named DNA, and the string of nucleotides making up the DNA should be passed as the sole argument 
# to the DNA initializer.
# Define a validate method that checks whether the letters given to the DNA object are valid nucleotides
# Add a __str__ method to your class, print in all uppercase


# **** TASK 2 ****
# Add a method convert_to_rna to the DNA class that returns a string containing the same data as the DNA sequence, 
# but with the Ts replaced by Us.
# Add a method reverse_complement to your DNA class. It should return a string that contains the reverse complement of the DNA sequence.


# **** TASK 3 ****
# Create an RNA class to model a strand of RNA.
# Similar to your DNA class, your RNA initializer should take and store a string representing the RNA bases. Your RNA class should
# have a validate method that returns whether all letters in the RNA sequence are valid nucleotides. Remember that this 
# is similar to DNA, but that the valid nucleotides in RNA are A, C, G, and U.
# Your class should also have a __str__ method that returns the uppercase RNA string (again similar to your DNA class).


# **** TASK 4 ****
# Translate an RNA sequence into a protein.

# ***Translation consists of scanning through the nucleotides, translating sets of 3 consecutive RNA bases into their corresponding codons. First, translation does not actually begin until a special start codon is encountered -- the codon M ('methionine'), in particular. Translation continues until a stop codon is encountered -- the stop codons are UAG ('amber'), UAA ('ochre'), and UGA ('opal' or 'umber')***
# Write a method to_protein inside your RNA class that will translate the RNA string into a protein, returning a string 
# corresponding to the protein. The method will take one argument -- the reading frame -- which should be either the value 0, 1, or 2.


# **** TASK 4 ****
# You will read in a file containing DNA sequence data, construct the corresponding DNA objects, convert those to RNA objects, 
# and finally convert the RNA objects to proteins.
# Your answer must include your read_dna_file function as well as your complete DNA and RNA classes.


