import global_align as g_a


keys = ['A', 'C', 'T', 'G', '-']
delta = {}
for i in range(len(keys)):
    delta[keys[i]] = {k: v for (k, v) in zip(keys, [1 if keys[i] == keys[j] else -1 for j in range(len(keys))])}


#Change the value to test
#---------------------------------#
A = "TAGATA"
B = "GTAGGCTTAAGGTTA"

# A = "T" * 1000
# B = "G" * 1000
#---------------------------------#

print('------------------Global Alignment-----------------')
score, align, space = g_a.global_align(A, B, delta)
aligned_A, aligned_B = align.split('\n')
print("Aligned A:", aligned_A)
print("Aligned B:", aligned_B)
print("Aligned Score:", score)
print("Space used:", space)
print('---------Global Alignment with Hirschberg----------')
aligned_A, aligned_B, score, space = g_a.hirschberg(A, B, delta)
print("Aligned A:", aligned_A)
print("Aligned B:", aligned_B)
print("Aligned Score:", score)
print("Space used:", space)
