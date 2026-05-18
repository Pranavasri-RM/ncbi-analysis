import os
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt


# Data Fetching

def fetch_sequence(accession_id, email, rettype="fasta"):
    """
    Retrieve a sequence record from the NCBI nucleotide database.

    Args:
        accession_id (str): NCBI accession ID (e.g. 'NM_000546').
        email (str): Email address required by NCBI Entrez API.
        rettype (str, optional): Return format - 'fasta' or 'gb' (GenBank).
            Defaults to 'fasta'.

    Returns:
        SeqRecord: Biopython SeqRecord containing the sequence and metadata.
    """
    Entrez.email = email
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession_id,
        rettype=rettype,
        retmode="text"
    )
    record = SeqIO.read(handle, rettype)
    handle.close()
    return record


# Sequence Statistics

def calculate_gc_content(record):
    """
    Calculate the GC content of a sequence as a percentage.

    Args:
        record (SeqRecord): Biopython SeqRecord object.

    Returns:
        float: GC content percentage, rounded to 2 decimal places.
    """
    seq = str(record.seq).upper()
    gc_count = seq.count("G") + seq.count("C")
    return round((gc_count / len(seq)) * 100, 2)



# ORF Detection

def find_orfs(record, min_length=100):
    """
    Identify Open Reading Frames (ORFs) in the forward strand of a sequence.

    Note: Nested/overlapping ORFs are included, as this is a naive finder.

    Args:
        record (SeqRecord): Biopython SeqRecord object.
        min_length (int, optional): Minimum ORF length in nucleotides.
            Defaults to 100.

    Returns:
        list[dict]: List of ORFs, each represented as a dictionary with keys:
            - 'start'    (int): Start position in the sequence.
            - 'end'      (int): End position (inclusive of stop codon).
            - 'length'   (int): Length in nucleotides.
            - 'sequence' (Seq): The ORF nucleotide sequence.
    """
    seq = record.seq.upper()
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    for i in range(len(seq) - 2):
        if seq[i:i + 3] == "ATG":
            for j in range(i, len(seq) - 2, 3):
                codon = str(seq[j:j + 3])
                if codon in stop_codons:
                    start, end = i, j + 3
                    length = end - start
                    if length >= min_length:
                        orfs.append({
                            "start": start,
                            "end": end,
                            "length": length,
                            "sequence": seq[start:end],
                        })
                    break  # stop at first in-frame stop codon

    return orfs



# Visualization

def plot_orf_distribution(orfs, output_path="results/orf_distribution.png"):
    """
    Plot a histogram of ORF lengths with the longest ORF marked.

    Args:
        orfs (list[dict]): ORF list returned by find_orfs().
        output_path (str, optional): File path for the saved figure.
            Defaults to 'results/orf_distribution.png'.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    lengths = [orf["length"] for orf in orfs]
    longest = max(lengths)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(lengths, bins=30, color="steelblue", edgecolor="black", alpha=0.85)
    ax.axvline(
        x=longest,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Longest ORF ({longest} nt)"
    )
    ax.set_xlabel("ORF Length (nucleotides)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("ORF Length Distribution — TP53 (NM_000546)", fontsize=14)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Plot saved to {output_path}")


# Main

if __name__ == "__main__":
    EMAIL = "your@email.com"
    ACCESSION = "NM_000546"

    print(f"Fetching sequence: {ACCESSION}")
    record = fetch_sequence(ACCESSION, EMAIL)
    print(f"ID          : {record.id}")
    print(f"Description : {record.description}")
    print(f"Length      : {len(record.seq)} nt")

    gc = calculate_gc_content(record)
    print(f"GC Content  : {gc}%")

    orfs = find_orfs(record)
    longest = max(orfs, key=lambda x: x["length"])
    print(f"\nORF Analysis:")
    print(f"Total ORFs found : {len(orfs)}")
    print(f"Longest ORF      : {longest['length']} nt (positions {longest['start']}-{longest['end']})")

    plot_orf_distribution(orfs)