#!/usr/bin/env python3
import subprocess
import os

# Annealing, Identifying Amplicon Pairs, and Extracting Amplicons, combined into a PCR function
def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    # Construct the BLASTN command
    blastn_cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", primer_file,
        "-subject", assembly_file,
        "-outfmt", "6 std qlen",
        "-evalue", "0.1"
    ]

    try:
        # Run the command and capture the output
        result = subprocess.run(blastn_cmd, capture_output=True, text=True, check=True)

        output_lines = result.stdout.strip().split('\n')
        # Create a list of matching hits that have a percent identity > 80%
        sorted_good_hits = []
        for line in output_lines:
            if line:
                fields = line.split('\t')
                percent_identity = float(fields[2])
                if percent_identity > 80.0:
                    sorted_good_hits.append(fields)
        # Sort the results by the sstart field (ninth field)
        sorted_good_hits.sort(key=lambda x: int(x[8]))
        hit_pairs = []
        i = 0
        n = len(sorted_good_hits)
        while i < n - 1:
            first_hit = sorted_good_hits[i]
            second_hit = sorted_good_hits[i + 1]

            # Assume sstart is in the ninth field
            sstart1 = int(first_hit[8])
            sstart2 = int(second_hit[8])

            # Check if the hits are within the max_amplicon_size
            if sstart2 - sstart1 < max_amplicon_size:
                # If they form a pair, add to the list
                hit_pairs.append((first_hit, second_hit))
                i += 2  # Move 2 indices because we already added a pair. Assume that each hit can only pair up with ONE other hit
            else:
                i += 1  # Move 1 index since no pair was found
        extracted_sequences = []
        for pair in hit_pairs:
            sseqid = pair[0][1]
            start = pair[0][9]  # Assuming ssend is in the tenth field
            end = str(int(pair[1][9]) - 1)  # Adjust for BED format
            bed_entry = f"{sseqid}\t{start}\t{end}\n"
            extracted_sequences.append(bed_entry)
        # Create a temporary BED file using mktemp
        bed_file = os.popen('mktemp /tmp/extracted_sequences.XXXXXX.bed').read().strip()

        # Write the BED entries to the temporary file
        with open(bed_file, 'w') as f:
            f.writelines(extracted_sequences)
        # Call seqtk to extract sequences
        cmd = ["seqtk", "subseq", assembly_file, bed_file]
        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.stdout

    except subprocess.CalledProcessError as e:
        print("Error running blastn:", e)
        return []


