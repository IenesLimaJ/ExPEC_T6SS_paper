#####################################################################
# December 6, 2024 - Julia Ienes-Lima

import subprocess
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import os

# Define the input files
query_file = "/path/to/query/sequence/T6SStemplate.fasta"
query_name = 'T6SS1/T6SS2/T6SS3' #change according to the query file name
pathotype = 'APEC/NMEC/UPEC/SEPEC' #change according to the pathotype that is being screening
subject_file = "/path/to/subject/sequence/genome.fasta"
output_file = "blast_output.xml"  # BLAST output in XML format
filtered_fasta = "filtered_hits.fasta"  # File to store the filtered sequences


# Run BLASTn with XML output format
def run_blast(query, subject, output): 
    blastn_command = [
        "blastn",
        "-query", query,
        "-subject", subject,
        "-outfmt", "5",  # XML format
        "-out", output
    ]
    subprocess.run(blastn_command)
    print(f"BLAST XML output saved to {output}")


# Parse BLAST XML output and filter for hits above a minimum length of 500bp.
def parse_blast_output(xml_file, min_length=500): 
    genome_sequences = defaultdict(str) # Use defaultdict to store sequences by genome ID
    with open(xml_file) as blast_output: # open the xml_file and ensure it's properly closed
        blast_records = NCBIXML.parse(blast_output)
        for record in blast_records:
            add_filtered_sequences(record, genome_sequences, min_length)
    return genome_sequences


# Add sequences from a BLAST record to genome_sequences if they exceed min_length (500 bp)
def add_filtered_sequences(blast_record, genome_sequences, min_length):
    for alignment in blast_record.alignments:
        genome_id = alignment.hit_id + pathotype + query_name 
        filtered_sequence = concatenate_long_sequences(alignment.hsps, min_length)
        if filtered_sequence:
            genome_sequences[genome_id] += filtered_sequence # associate the genome ID with the filtered sequences (hits in subject file)


# Concatenate HSP (High-scoring Segment Pairs) sequences that meet the minimum length requirement
# HSP = hit between the query and the subject sequence
def concatenate_long_sequences(hsps, min_length):
    return ''.join(hsp.sbjct for hsp in hsps if hsp.align_length > min_length) # merge the hits > 500bp from subject sequence into one


# Save concatenated sequences to a FASTA file
def save_to_fasta(genome_sequences, output_fasta):
    concatenated_sequences = [
        SeqRecord(Seq(sequence), id=genome_id, description="Concatenated sequence")
        for genome_id, sequence in genome_sequences.items()
    ]
    SeqIO.write(concatenated_sequences, output_fasta, "fasta")
    print(f"Concatenated filtered sequences saved to {output_fasta}")

# Run BLASTn and generate XML output
run_blast(query_file, subject_file, output_file)

# Parse BLAST XML and save concatenated sequences
genome_sequences = parse_blast_output(output_file)
save_to_fasta(genome_sequences, filtered_fasta)

#####################################################################
