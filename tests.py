import numpy as np
from sphertree import Tree

def test_construction():
    arr = np.array([[1., 2., 3.], [4., 5., 6.]], dtype='f4')
    tree = Tree(arr)
    assert(tree[0] == (1., 2., 3., 0))
    assert(tree[1] == (4., 5., 6., 1))

if __name__ == "__main__":
    test_construction()
