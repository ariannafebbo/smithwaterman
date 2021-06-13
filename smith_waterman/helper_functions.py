import numpy as np

code2dir = {
  0: 'D',
  1: 'L',
  2: 'U',
  3: 'N',
}

def matrix_score(seq1, seq2, match_cost=1, mismatch_cost=-1, gap_cost=-2):
    S=np.zeros((len(seq1)+1, len(seq2)+1), dtype=int)
    D=np.full((len(seq1)+1, len(seq2)+1), 'N', dtype=str)
    
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if seq1[i-1] == seq2[j-1]:
                diagonal_score = S[i-1, j-1] + match_cost
            else:
                diagonal_score = S[i-1, j-1] + mismatch_cost
            left_score = S[i, j-1] + gap_cost
            up_score = S[i-1, j] + gap_cost
            D[i, j] = code2dir[np.argmax([diagonal_score, left_score, up_score, 0])]
            S[i, j] = np.max([diagonal_score, left_score, up_score, 0])
    return S, D

def traceback(score_matrix, direction_matrix):
  i,j = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
  max_score = score_matrix[i,j]
  path = []

  while direction_matrix[i,j] != 'N':
    path.append([i,j])

    if direction_matrix[i,j] == 'D':
      i = i - 1
      j = j - 1
    elif direction_matrix[i,j] == 'L':
      j = j - 1
    elif direction_matrix[i,j] == 'U':
      i = i - 1

  return path, max_score

def path2alignment(seq1, seq2, path):
  align1 = ""
  align2 = ""

  path.reverse()

  seq1_index = path[0][0]
  seq2_index = path[0][1]
  align1 += seq1[seq1_index - 1]
  align2 += seq2[seq2_index - 1]

  for ind in path[1:]:
    if (ind[0] == seq1_index):
      align1 += '-'
    else:
      seq1_index = ind[0]
      align1 += seq1[seq1_index-1]
        
    if (ind[1] == seq2_index):
      align2+='-'
    else:
      seq2_index = ind[1]
      align2+=seq2[seq2_index-1]

  return align1, align2