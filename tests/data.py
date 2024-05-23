

VALID_SCHEMA = {
    'A': 0.31, 'L': 1.70, 'R': -1.01, 'K': -0.99, 'N': -0.60,
    'M': 1.23, 'D': -0.77, 'F': 1.79, 'C': 1.54, 'P': 0.72,
    'Q': -0.22, 'S': -0.04, 'E': -0.64, 'T': 0.26, 'G': 0.0,
    'W': 2.25, 'H': 0.13, 'Y': 0.96, 'I': 1.80, 'V': 1.22
}
INVALID_SCHEMA = {
    'A': 0.31, 'L': 1.70, 'R': -1.01, 'K': -0.99, 'N': -0.60,
    'M': 1.23, 'D': -0.77, 'F': 1.79, 'C': 1.54, 'P': 0.72,
    'Q': -0.22, 'S': -0.04, 'E': -0.64, 'T': 0.26, 'G': 0.0,
    'W': 2.25, 'H': 0.13, 'Y': 0.96, 'I': 1.80
}

VALID_SEQUENCES = [
    "MFVFLVLLPLVSSQCVNL",
    "FKIYSKHTP",
    "QDVNCTEVPVAIHADQLTPTWRVY"
]
INVALID_SEQUENCES = [
    "XMFVFLVLLPLVSSQCVNL",
    "ZFKIYSKHTP",
    ".QDVNCTEVPVAIHADQLTPTWRVY"
]

