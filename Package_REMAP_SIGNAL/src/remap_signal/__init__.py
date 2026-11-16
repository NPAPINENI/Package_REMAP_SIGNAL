from .core import RemapSignal, SPECIES_ALIASES

__all__ = ["RemapSignal", "SPECIES_ALIASES"]
__version__ = "0.1.0"

from .pipeline import build_datasets, process_tf, process_file, make_splits  # noqa: E402
__all__ += ["build_datasets", "process_tf", "process_file", "make_splits"]

from .multi import Remap, load, resolve_bigwig  # noqa: E402
__all__ += ["Remap", "load", "resolve_bigwig"]

from .download import download  # noqa: E402
__all__ += ["download"]

from .api import get_signal, add_signal  # noqa: E402
__all__ += ["get_signal", "add_signal"]
