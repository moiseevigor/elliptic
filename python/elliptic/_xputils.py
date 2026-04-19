"""Array-namespace detection helper shared across all modules."""
from __future__ import annotations

from array_api_compat import array_namespace, is_array_api_obj
import numpy as np


def get_xp(*args):
    """Return the array namespace for *args*, defaulting to numpy for plain scalars."""
    api_objs = [a for a in args if is_array_api_obj(a)]
    if api_objs:
        return array_namespace(*api_objs)
    return np
