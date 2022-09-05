import urllib.request as ur
from Bio import SeqIO
import numpy as np
from tqdm import tqdm
import pokemons as pnames
import xml.etree.ElementTree as ET


def generate_fasta_file_with_names(path="./hemoglobin.fasta", N=10):
    """
        Use file and instructions provided at 
        https://www.mimuw.edu.pl/~lukaskoz/teaching/sad2/lab6/readme.html
    """
    records = list(SeqIO.parse(path, "fasta"))
    # Only pick ones that are in the 137-143 length range
    wanted_records = filter(lambda x: len(x.seq) >= 139 and len(x.seq) <= 141, records)
    wanted_records = list(wanted_records)[:N]
    for record in tqdm(wanted_records, "Fetching information"):
        full_name, short_name, group = fetch_identifier_info(record)
        record.description = "|" + full_name + "|" + short_name + "|" + group
    # Write relevant information to a new file
    SeqIO.write(wanted_records, "hemoglobin_with_names.fasta", "fasta")


def fetch_identifier_info(record):
    base_url = "https://www.uniprot.org/uniprot/"
    identifier = record.id
    url = base_url + identifier + ".xml"
    with ur.urlopen(url) as u:
        # Parse XML content
        tree = ET.ElementTree(ET.fromstring(u.read()))
    # Extract name and shortSpeciesName
    namespace = "{http://uniprot.org/uniprot}"
    root = tree.getroot()
    wanted = root.find(namespace + "entry").find(namespace + "organism")
    for child in wanted:
        if child.get('type') == 'scientific':
            scientific_name = child.text
    name = scientific_name
    shortSpeciesName = name.split(" ")[0][0] + name.split(" ")[1][:3]
    # Get lineage as well
    group_info = wanted.find(namespace + "lineage").findall(namespace + "taxon")
    group = group_info[6].text
    return name, shortSpeciesName, group


def get_sequences_for_tree():
    records = list(SeqIO.parse("hemoglobin_with_names.fasta", "fasta"))
    processed_records = []
    for record in records:
        sequence = record.seq
        id = record.id
        fullname, name, group = record.description.split("|")[1:]
        processed_records.append(((id, fullname, name, group), str(sequence)))
    return processed_records


def make_sequences_for_ancestry(length = 15, n_seq = 30):
    sequences = []

    for i in range(n_seq):
        # Generate random nucleotide sequence
        sequence = np.random.choice(4, size=length)
        sequence ="".join(["ATCG"[x] for x in sequence])
        sequences.append(sequence)

    names = pnames.genarr(len(sequences))
    # Save these to a file
    with open("ancestry_info.txt", "w") as f:
        for i in range(len(sequences)):
            f.write(sequences[i] + "," + names[i] + "\n")


def get_sequences_for_ancestry():
    sequences, names = [], []
    with open("ancestry_info.txt", "r") as f:
        for line in f:
            sequence, name = line.strip("\n").split(",")
            sequences.append(sequence)
            names.append(name)
    return sequences, names


def convert_to_amino(seq):
    codons = {
        "UUU": "F",
        "CUU": "L",
        "AUU": "I",
        "GUU": "V",
        "UUC": "F",
        "CUC": "L",
        "AUC": "I",
        "GUC": "V",
        "UUA": "L",
        "CUA": "L",
        "AUA": "I",
        "GUA": "V",
        "UUG": "L",
        "CUG": "L",
        "AUG": "M",
        "GUG": "V",
        "UCU": "S",
        "CCU": "P",
        "ACU": "T",
        "GCU": "A",
        "UCC": "S",
        "CCC": "P",
        "ACC": "T",
        "GCC": "A",
        "UCA": "S",
        "CCA": "P",
        "ACA": "T",
        "GCA": "A",
        "UCG": "S",
        "CCG": "P",
        "ACG": "T",
        "GCG": "A",
        "UAU": "Y",
        "CAU": "H",
        "AAU": "N",
        "GAU": "D",
        "UAC": "Y",
        "CAC": "H",
        "AAC": "N",
        "GAC": "D",
        "UAA": "_",
        "CAA": "Q",
        "AAA": "K",
        "GAA": "E",
        "UAG": "_",
        "CAG": "Q",
        "AAG": "K",
        "GAG": "E",
        "UGU": "C",
        "CGU": "R",
        "AGU": "S",
        "GGU": "G",
        "UGC": "C",
        "CGC": "R",
        "AGC": "S",
        "GGC": "G",
        "UGA": "_",
        "CGA": "R",
        "AGA": "R",
        "GGA": "G",
        "UGG": "W",
        "CGG": "R",
        "AGG": "R",
        "GGG": "G"
    }
    seq = seq.replace("T", "U")
    amino = []
    for i in range(0, len(seq), 3):
        amino.append(codons[seq[i:i+3]])
    return "".join(amino)

def load_oca2_sequences():
    # Human ends with "TAAAA", of which "TAA" is a stop codon
    # Alter sequence to remove last two nucleotides (for compatibility with amino-acid question)
    human_oca2 = str(
        list(SeqIO.parse("human_oca2.fasta", "fasta"))[0].seq)
    human_oca2 = human_oca2[:-2]
    # Mouse ends with "TA" and is not divislbe by 3
    # Good chance the last character is A/G, which would make it a stop codon
    mouse_oca2 = str(
        list(SeqIO.parse("mouse_oca2.fasta", "fasta"))[0].seq)
    mouse_oca2+= "A"
    # Crop sequences (too long currently)
    human_oca2 = human_oca2[:219]
    mouse_oca2 = mouse_oca2[:207]

    return human_oca2, mouse_oca2
    
    
def get_hemoglobin_sequences():
    polar_bear = str(
        list(SeqIO.parse("polar_bear_hemoglobin.fasta", "fasta"))[0].seq)
    black_bear = str(
        list(SeqIO.parse("black_bear_hemoglobin.fasta", "fasta"))[0].seq)
    human = str(
        list(SeqIO.parse("human_hemoglobin.fasta", "fasta"))[0].seq)
    chimp = str(
        list(SeqIO.parse("chimp_hemoglobin.fasta", "fasta"))[0].seq)
    return polar_bear,  black_bear, human, chimp
