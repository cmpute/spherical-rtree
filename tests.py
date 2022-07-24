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

if __name__ == "__main__":
    test_construction()
    test_initialize()
