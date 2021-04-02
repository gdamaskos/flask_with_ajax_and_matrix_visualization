import os
from contextlib import contextmanager


@contextmanager
def env_var(key, value):
    # Save old value, if it exists
    old = None
    if key in os.environ:
        old = os.environ[key]
    # Set (new) env var
    os.environ[key] = value
    yield
    # Restore or delete env var
    if old:
        os.environ[key] = old
    else:
        del os.environ[key]
