from typing import List, Tuple

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

Predicates = Literal['cone', 'cone_frustum', 'ball']

class Tree:
    """
    Spherical RTree instance
    """
    def __init__(self, points: List[List[float]], values: List[int]) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Tuple[float, float, float, int]: ...
    def initialize(self, origin: Tuple[float, float, float]) -> None: ...
    def dispose(self) -> None: ...
    def is_initialized(self) -> bool: ...
    def origin(self) -> Tuple[float, float, float]: ...
    def query(self, point: Tuple[float, float, float], pred: Predicates, pred_args: int) -> List[int]: ...