import sys
import numpy as np
from numpy import linalg as LA

node_count = 0
alpha = 0
beta = 0

# Wheel graph is symmetric so it suffices to iterate over lower diagonal
def wheelgraph(n):
	mat = np.zeros([n,n])
	for i in range(n):
		for j in range(i):
			if (i == j):
				continue          
			elif i == 0 or j == 0:
				mat[i,j] = 1
				mat[j,i] = 1
			else:
				left_neighbor = i - 1
				right_neighbor = i + 1
				
				if left_neighbor < 0:
					left_neighbor = n-1
				if right_neighbor >= n:
					right_neighbor = (right_neighbor + 1) % n

				print (left_neighbor, i, right_neighbor)
				mat[i, left_neighbor] = 1
				mat[i, right_neighbor] = 1
				mat[left_neighbor, i] = 1
				mat[right_neighbor, i] = 1
	return mat

def katz_centrality(graph, alph, bet):
	print graph
	column_ones = np.ones(len(graph)).reshape(len(graph), 1)
	inverse_interm = inverse_intermediate(graph, alph)
	inverse_interm_beta = bet * inverse_interm
	return np.dot(inverse_interm_beta, column_ones)

def inverse_intermediate(graph, alph):
	identity = np.identity(len(graph))
	interm_mat = identity - (alph * np.transpose(graph))
	return LA.inv(interm_mat)

def get_eigenvalues(graph):
	return LA.eig(graph)

# Defaults to 6 node graph, with alpha = 0.2, and beta = 1.0
if len(sys.argv) == 1:
	node_count = 6
	alpha = 0.2
	beta = 1.0
elif len(sys.argv) == 2:
	node_count = int(sys.argv[1])
	alpha = 1.0/(node_count - 1)
	beta = 1.0
elif len(sys.argv) == 3:
	node_count = int(sys.argv[1])
	alpha = float(sys.argv[2])
	beta = 1.0
elif len(sys.argv) == 4:
	node_count = int(sys.argv[1])
	alpha = float(sys.argv[2])
	beta = float(sys.argv[3])

wheel = wheelgraph(node_count)
print katz_centrality(wheel, alpha, beta)