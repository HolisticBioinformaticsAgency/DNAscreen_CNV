from Bio import SeqIO

# ---- INPUT / OUTPUT ----
input_fasta  = "/fs04/vh83/jason/ont/barcodes_RC_F.fa"
output_fasta = "/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/demultiplex_pb/workflow_HiFiTargetEnrichment/barcodes/custom_lima_barcodes_v2.fasta"

# ---- ADAPTERS (from your Dorado TOML) ----
MASK1_FRONT = "AATGATACGGCGACCACCGAGATCTACAC"
MASK1_REAR  = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"

MASK2_FRONT = "CAAGCAGAAGACGGCATACGAGAT"
MASK2_REAR  = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

# ---- PROCESS ----
with open(output_fasta, "w") as out:
    for record in SeqIO.parse(input_fasta, "fasta"):
        name = record.id
        seq  = str(record.seq).upper()

        if name.startswith("BCF"):
            bc_id = name.replace("BCF", "bc", 1)
            header = f">{bc_id}_r"
            full_seq = MASK2_FRONT + seq + MASK2_REAR

        elif name.startswith("BCR"):
            bc_id = name.replace("BCR", "bc", 1)
            header = f">{bc_id}_f"
            full_seq = MASK1_FRONT + seq + MASK1_REAR

        else:
            raise ValueError(f"Unrecognized barcode name: {name}")

        out.write(header + "\n")
        out.write(full_seq + "\n")

print(f"Wrote lima-compatible barcode file: {output_fasta}")
