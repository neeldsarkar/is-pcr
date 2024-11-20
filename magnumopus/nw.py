#!/usr/bin/env python3

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    # Initialize the scoring matrix
    len_a = len(seq_a) + 1
    len_b = len(seq_b) + 1
    score_matrix = [[0] * len_b for _ in range(len_a)]

    # Fill in the first row and first column of the scoring matrix
    for i in range(len_a):
        score_matrix[i][0] = gap * i
    for j in range(len_b):
        score_matrix[0][j] = gap * j

    # Fill in the rest of the scoring matrix
    for i in range(1, len_a):
        for j in range(1, len_b):
            match_score = score_matrix[i-1][j-1] + (match if seq_a[i-1] == seq_b[j-1] else mismatch)
            delete_score = score_matrix[i-1][j] + gap
            insert_score = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(match_score, delete_score, insert_score)

    # Backtrack to find the optimal alignment
    aligned_a, aligned_b = [], []
    i, j = len_a - 1, len_b - 1

    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + (match if seq_a[i-1] == seq_b[j-1] else mismatch):
            aligned_a.append(seq_a[i-1])
            aligned_b.append(seq_b[j-1])
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i-1][j] + gap:
            aligned_a.append(seq_a[i-1])
            aligned_b.append('-')
            i -= 1
        else:
            aligned_a.append('-')
            aligned_b.append(seq_b[j-1])
            j -= 1

    # Reverse the aligned sequences since they were constructed backwards
    aligned_a.reverse()
    aligned_b.reverse()

    # Calculate the final score
    final_score = score_matrix[len_a - 1][len_b - 1]

    return (''.join(aligned_a), ''.join(aligned_b)), final_score
