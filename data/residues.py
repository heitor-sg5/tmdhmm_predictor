# Residue hydrophobicity (KD) scores
KD_SCORES = {
    "A":  1.8, # Alanine
    "R": -4.5, # Arginine
    "N": -3.5, # Asparagine
    "D": -3.5, # Aspartate
    "C":  2.5, # Cysteine
    "Q": -3.5, # Glutamine
    "E": -3.5, # Glutamate  
    "G": -0.4, # Glycine
    "H": -3.2, # Histidine 
    "I":  4.5, # Isoleucine 
    "L":  3.8, # Leucine
    "K": -3.9, # Lysine
    "M":  1.9, # Methionine
    "F":  2.8, # Phenylalanine 
    "P": -1.6, # Proline
    "S": -0.8, # Serine 
    "T": -0.7, # Threonine
    "W": -0.9, # Tryptophan 
    "Y": -1.3, # Tyrosine  
    "V":  4.2, # Valine 
}

# Residue classes
RESIDUE_CLASS = {
    # Strongly hydrophobic (TM)
    "I": "H", "L": "H", "V": "H", 
    "F": "H", "A": "H", "M": "H",
    # Aromatic (caps)
    "W": "A", "Y": "A",
    # Polar uncharged (aqueous)
    "S": "P", "T": "P", "N": "P", 
    "Q": "P", "C": "P", "G": "P",
    # Charged (aqueous)
    "R": "Q", "K": "Q", 
    "D": "Q", "E": "Q", "H": "Q",
    # Proline (helix breaker)
    "P": "P",
}

# Residues enriched on cytosolic side (positive-inside rule)
POSITIVE_RESIDUES = frozenset("RK")

# Aromatic residues that anchor TM helices at the membrane interface
AROMATIC_RESIDUES = frozenset("WYF")

# All valid single-letter amino acid codes
VALID_AA = frozenset(KD_SCORES.keys())

# Full names
AA_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartate",
    "C": "Cysteine", "Q": "Glutamine", "E": "Glutamate", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
}