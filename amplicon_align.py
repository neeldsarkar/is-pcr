#!/usr/bin/env python3
import magnumopus
import argparse

parser = argparse.ArgumentParser(
    description="Perform in-silico PCR on two assemblies and align the amplicons.",
    formatter_class=argparse.RawTextHelpFormatter
)

# Add arguments
parser.add_argument('-1', '--assembly1', type=str, required=True, help='Path to the first assembly file.')
parser.add_argument('-2', '--assembly2', type=str, required=True, help='Path to the second assembly file.')
parser.add_argument('-p', '--primers', type=str, required=True, help='Path to the primer file.')
parser.add_argument('-m', '--max_amplicon_size', type=int, required=True, help='Maximum amplicon size for isPCR.')
parser.add_argument('--match', type=int, default=1, help='Match score to use in alignment.')
parser.add_argument('--mismatch', type=int, default=-1, help='Mismatch penalty to use in alignment.')
parser.add_argument('--gap', type=int, default=-1, help='Gap penalty to use in alignment.')

args = parser.parse_args()

first_match = magnumopus.ispcr(args.primers, args.assembly1, args.max_amplicon_size)
second_match = magnumopus.ispcr(args.primers, args.assembly2, args.max_amplicon_size)

seq1 = ''.join(line.strip() for line in first_match.splitlines() if not line.startswith('>'))
seq2 = ''.join(line.strip() for line in second_match.splitlines() if not line.startswith('>'))
alns, score = magnumopus.needleman_wunsch(seq1, seq2, 1, -1, -1)
for aln in alns:
    print(aln)
print(score)

