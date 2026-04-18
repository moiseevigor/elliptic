"""Pytest configuration: backend fixtures."""
import pytest
import numpy as np

BACKENDS = ["numpy"]

try:
    import torch
    BACKENDS.append("torch")
except ImportError:
    pass

try:
    import jax.numpy  # noqa: F401
    BACKENDS.append("jax")
except ImportError:
    pass


@pytest.fixture(params=BACKENDS)
def xp(request):
    name = request.param
    if name == "numpy":
        return np
    if name == "torch":
        import torch
        return torch
    if name == "jax":
        import jax.numpy as jnp
        return jnp
    raise ValueError(name)


def to_array(xp, data):
    """Convert a numpy array to the given backend."""
    if xp is np:
        return np.asarray(data)
    name = type(xp).__name__
    if "torch" in str(xp):
        import torch
        return torch.tensor(data, dtype=torch.float64)
    import jax.numpy as jnp
    return jnp.array(data)
