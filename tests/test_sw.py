
from smith_waterman import sw
import numpy as np
import random 
import string

def test_find_alignment_stress():
  letters = string.ascii_lowercase
  seq_len = 500
  s1 = ''.join(random.choice(letters) for i in range(seq_len))
  align1, align2, s = sw.find_alignment(s1, s1, match_cost=3, mismatch_cost=-3, gap_cost=-2)
  assert align1 == align2

def test_find_alignment():
  seq1 = "GGTTGACTA"
  seq2 = "TGTTACGG"
  expected_align1 = "GTTGAC"
  expected_align2 = "GTT-AC"
  align1, align2, max_score = sw.find_alignment(seq1, seq2, match_cost=3, mismatch_cost=-3, gap_cost=-2)
  assert align1 == expected_align1, "Wrong alignment 1"
  assert align2 == expected_align2, "Wrong alignment 2"
  assert max_score == 13, "Wrong alignment score"
