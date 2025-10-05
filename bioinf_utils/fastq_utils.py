from typing import Union, Tuple

def normalize_bounds(bounds: Union[int, float, Tuple[float, float]]) -> Tuple[float, float]:
    """
    Normalize bounds to a (min, max) tuple.

    Arguments:
    bounds: single number (interpreted as upper bound) or pair (min, max)

    Returns tuple of two floats (min, max).
    Raises ValueError if input is invalid.
    """
    if isinstance(bounds, (int, float)):
        return (0.0, float(bounds))
    elif isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        return (float(bounds[0]), float(bounds[1]))
    else:
        raise ValueError("Bounds must be a number or a pair of numbers")