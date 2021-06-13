
from smith_waterman import sw
import numpy as np
import random 
import string

def test_smith_waterman_stress():
  letters = string.ascii_lowercase
  result_str = ''.join(random.choice(letters) for i in range(500))
  align1, align2, s = sw.smithwaterman(result_str, result_str, match_cost=3, mismatch_cost=-3, gap_cost=-2)
  assert align1 == align2

def test_smith_waterman():
  seq1 = "GGTTGACTA"
  seq2 = "TGTTACGG"
  expected_align1 = "GTTGAC"
  expected_align2 = "GTT-AC"
  align1, align2, max_score = sw.smithwaterman(seq1, seq2, match_cost=3, mismatch_cost=-3, gap_cost=-2)
  assert align1 == expected_align1, "Wrong alignment 1"
  assert align2 == expected_align2, "Wrong alignment 2"
  assert max_score == 13, "Wrong alignment score"

def test_path2align():
  seq1 = "ACCG"
  seq2 = "CACC"
  path = [[4,4], [4,3], [3,3], [2,2], [2,1]]
  expected_align1 = "C-CG-"
  expected_align2 = "CAC-C"
  align1, align2 = sw.path2alignment(seq1, seq2, path)
  assert align1 == expected_align1, "Alignment 1 is wrong"
  assert align2 == expected_align2, "Alignment 2 is wrong"

def test_traceback_seq():
  score_matrix = np.array([
    [0, 0, 0, 0, 0],
    [0, 1, 3, 1, 1],
    [0, 0, 3, 1, 3],
    [0, 3, 9, 6, 10],
    [0, 2, 1, 3, 1]
    ])
  direction_matrix = np.array([
    ['N', 'N', 'N', 'N', 'N'],
    ['N', 'D', 'D', 'L', 'U'],
    ['N', 'N', 'D', 'L', 'D'],
    ['N', 'D', 'L', 'D', 'L'],
    ['N', 'D', 'L', 'L', 'L']
    ])
  expected_result = np.array([[3,4], [3,3], [2,2], [1,1]])
  path, max_score = sw.traceback(score_matrix, direction_matrix)
  assert (path == expected_result).all(), "Traceback path is wrong"
  assert max_score == 10, "Max score is wrong"

def test_long_seq():
    seq1 = "GGTTGACTA"
    seq2 = "TGTTACGG"
    expected_score_matrix = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 3, 1, 0, 0, 0, 3, 3],
                               [0, 0, 3, 1, 0, 0, 0, 3, 6],
                               [0, 3, 1, 6, 4, 2, 0, 1, 4],
                               [0, 3, 1, 4, 9, 7, 5, 3, 2],
                               [0, 1, 6, 4, 7, 6, 4, 8, 6],
                               [0, 0, 4, 3, 5, 10, 8, 6, 5],
                               [0, 0, 2, 1, 3, 8, 13, 11, 9],
                               [0, 3, 1, 5, 4, 6, 11, 10, 8],
                               [0, 1, 0, 3, 2, 7, 9, 8, 7]])
    expected_direction_matrix = np.array([['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'],
                                          ['N', 'N', 'D', 'L', 'N', 'N', 'N', 'D', 'D'],
                                          ['N', 'N', 'D', 'L', 'N', 'N', 'N', 'D', 'D'],
                                          ['N', 'D', 'L', 'D', 'D', 'L', 'L', 'U', 'U'],
                                          ['N', 'D', 'L', 'D', 'D', 'L', 'L', 'L', 'U'],
                                          ['N', 'U', 'D', 'L', 'U', 'D', 'D', 'D', 'D'],
                                          ['N', 'N', 'U', 'D', 'U', 'D', 'L', 'L', 'D'],
                                          ['N', 'N', 'U', 'D', 'U', 'U', 'D', 'L', 'L'],
                                          ['N', 'D', 'L', 'D', 'D', 'U', 'U', 'D', 'D'],
                                          ['N', 'U', 'D', 'U', 'D', 'D', 'U', 'D', 'D']])
    result_score_matrix, result_direction_matrix = sw.matrix_score(seq1, seq2, match_cost=3, mismatch_cost=-3, gap_cost=-2)
    assert (result_score_matrix == expected_score_matrix).all(), "Score matrix is not correct"
    assert (result_direction_matrix == expected_direction_matrix).all(), "Direction matrix is not correct"

def test_empty_seq():
  seq1 = "GGTTGACTA"
  seq2 = ""
  expected_score_matrix = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0]])
  expected_direction_matrix = np.array([['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N']])
  result_score_matrix, result_direction_matrix = sw.matrix_score(seq1, seq2)
  assert (result_score_matrix == expected_score_matrix).all(), "Score matrix is not correct"
  assert (result_direction_matrix == expected_direction_matrix).all(), "Direction matrix is not correct"