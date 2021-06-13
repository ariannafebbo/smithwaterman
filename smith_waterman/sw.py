from .helper_functions import compute_matrix_score, traceback, path2alignment

def find_alignment(seq1, seq2, match_cost=1, mismatch_cost=-1, gap_cost=-2):
  """
  Returns the optimal local alignment for the two input sequences and the alignment score,
  given the sequences I want to align and the cost of a match, of a mismatch
  and of a gap.

  Parameters
  ----------
  seq1 : str  
    first sequence to align
  seq2 : str
    second sequence to align
  match_cost : int
    score given when two bases match
  mismatch_cost : int
    score given when two bases do not match
  gap_cost : int
    score given when a base is aligned with a gap

  Returns
  -------
  align1, align2, max_score : tuple
    align1 is a string containing the bases alignment of seq1;
    align2 is a string containing the bases alignment of seq2;
    max_score is the alignment score
  
  """
  S,D = compute_matrix_score(seq1, seq2, match_cost, mismatch_cost, gap_cost)
  path, max_score = traceback(S,D)
  align1, align2 = path2alignment(seq1, seq2, path)
  return align1, align2, max_score
