from .helper_functions import matrix_score, traceback, path2alignment

def find_alignment(seq1, seq2, match_cost=1, mismatch_cost=-1, gap_cost=-2):
  S,D = matrix_score(seq1, seq2, match_cost, mismatch_cost, gap_cost)
  path, max_score = traceback(S,D)
  align1, align2 = path2alignment(seq1, seq2, path)
  return align1, align2, max_score
