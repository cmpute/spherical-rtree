import numpy as np
import time
from sphertree import Tree

def test_construction():
    arr = np.array([[1., 2., 3.], [4., 5., 6.]], dtype='f4')
    tree = Tree(arr)
    assert(tree[0] == (1., 2., 3., 0))
    assert(tree[1] == (4., 5., 6., 1))

def test_initialize():
    arr = np.random.rand(100, 3).astype('f4')
    tree = Tree(arr)
    tree.initialize((0, 0, 0))
    tree.initialize((1, 1, 1)) # override existing tree

def generate_box():
    N = 1000 # per face
    coords = np.random.rand(10, N) * 2 - 1
    right = np.array([np.ones(N), coords[0], coords[1]])
    left = np.array([-np.ones(N), coords[2], coords[3]])
    front = np.array([coords[4], -np.ones(N), coords[5]])
    back = np.array([coords[6], np.ones(N), coords[7]])
    up = np.array([coords[8], coords[9], np.ones(N)])
    sides = np.hstack([left, right, front, back, up])
    return sides.T.astype('f4')

def test_query():
    tree = Tree(generate_box())
    tree.initialize((0, 0, 0))
    assert(len(tree.query((1e-5, 1e-5, 2), "cone", 1)) > 0)
    assert(len(tree.query((2, 0, 0), "cone", 1)) > 0)
    assert(len(tree.query((1e-5, 1e-5, -2), "cone", 1)) == 0)

if __name__ == "__main__":
    test_construction()
    test_initialize()
    test_query()
