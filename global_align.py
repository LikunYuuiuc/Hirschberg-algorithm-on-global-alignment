import sys
UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)

def traceback_global(v, w, pointers):
    i,j = len(v), len(w)
    new_v = []
    new_w = []
    while True:
        di, dj = pointers[i][j]
        if (di,dj) == LEFT:
            new_v.append('-')
            new_w.append(w[j-1])
        elif (di,dj) == UP:
            new_v.append(v[i-1])
            new_w.append('-')
        elif (di,dj) == TOPLEFT:
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    return ''.join(new_v[::-1])+'\n'+''.join(new_w[::-1])



def global_align(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
    computed by traceback_global.

    :param: v
    :param: w
    :param: delta
    """
    # print("here: ", v)
    # print("here: ", w)
    # print("here: ", delta)
    M = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)]
    # print("here: ", M)
    # print("here: ", pointers)
    score, alignment = None, None
    # # YOUR CODE HERE
    # raise NotImplementedError()
    for i in range(1, len(v) + 1):
      M[i][0] = M[i-1][0] + delta[v[i-1]]['-']
      pointers[i][0] = UP

    for j in range(1, len(w) + 1):
      M[0][j] = M[0][j-1] + delta[w[j-1]]['-']
      pointers[0][j] = LEFT

    for i in range(1, len(v) + 1):
      for j in range(1, len(w) + 1):
        v_char = v[i-1]
        w_char = w[j-1]

        scoreul = M[i-1][j-1] + delta[v_char][w_char]
        scoreu = M[i-1][j] + delta[v_char]['-']
        scorel = M[i][j-1] + delta[w_char]['-']
        if scoreul >= scoreu and scoreul >= scorel:
          M[i][j] = scoreul
          pointers[i][j] = TOPLEFT
        elif scoreu > scoreul and scoreu > scorel:
          M[i][j] = scoreu
          pointers[i][j] = UP
        elif scorel > scoreul and scorel > scoreu:
          M[i][j] = scorel
          pointers[i][j] = LEFT


    score = M[len(v)][len(w)]

    space = sys.getsizeof(M) + sys.getsizeof(pointers)
    # print(score)
    # print(pointers)
    alignment = traceback_global(v,w, pointers)
    return score, alignment, space



def hirschberg(A, B, delta):
    """
    Hirschberg algorithm for sequence alignment.
    Args:
        A: First sequence (string or list).
        B: Second sequence (string or list).
    Returns:
        Aligned version of sequence A and B as a tuple.
    """
    max_space = float('-inf')
    if len(A) == 0:
        return ('-' * len(B), B, len(B), max_space)
    if len(B) == 0:
        return (A, '-' * len(A), len(A), max_space)
    if len(A) == 1 or len(B) == 1:
        score, align, space = global_align(A, B, delta)
        align = align.split('\n')
        max_space = max(max_space, space)
        return align[0], align[1], score, max_space

    mid = len(A) // 2

    # Compute score for left and right parts
    score_left, space = nw_score(A[:mid], B, delta)
    max_space = max(max_space, space)
    score_right, space = nw_score(A[mid:][::-1], B[::-1], delta)
    max_space = max(max_space, space)

    # Find the partition point

    max_score = float('-inf')
    partition = 0

    for j in range(len(B) + 1):
        current_score = score_left[j] + score_right[len(B) - j]
        if current_score > max_score:
            max_score = current_score
            partition = j

    # Recursively align the left and right parts
    left_A, left_B, left_score, max_space = hirschberg(A[:mid], B[:partition], delta)
    right_A, right_B, right_score, max_space = hirschberg(A[mid:], B[partition:], delta)



    return (left_A + right_A, left_B + right_B, left_score + right_score, max_space)



def nw_score(X, Y, delta):
    """
    Compute Needleman-Wunsch alignment scores for a single row.
    Args:
        X: Sequence X.
        Y: Sequence Y.
    Returns:
        A list containing alignment scores for the final row.
    """
    prev = []
    space = 0
    for i in range(len(Y) + 1):
        prev.append(0 - i)
    # prev = list(range(len(Y) + 1))
    for x in X:
        curr = [prev[0] - 1]
        for j in range(1, len(Y) + 1):
            match = prev[j - 1] + delta[x][Y[j - 1]]
            # here = prev[j]
            # here = delta[x]['-']
            delete = prev[j] + delta[x]['-']
            insert = curr[j - 1] + delta[Y[j-1]]['-']
            curr.append(max(match, delete, insert))
        prev = curr
        space = sys.getsizeof(curr)
    space += sys.getsizeof(prev)
    return prev, space

