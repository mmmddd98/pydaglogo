from importlib.metadata import PackageNotFoundError, version as _metadata_version

try:
    __version__ = _metadata_version("pydaglogo")
except PackageNotFoundError:
    __version__ = "unknown"

from .background import build_background_model, build_background_model_from_samples, load_background_samples
from .io import build_windows_from_flanks, build_windows_from_positions, format_sequence, prepare_proteome
from .models import BackgroundModel, DAUResult, DagPeptides, GroupingScheme, Proteome
from .plotting import plot_logo
from .statistics import test_dau, to_logo_matrix
from .workflow import make_daglogo, run_daglogo

__all__ = [
    "__version__",
    "BackgroundModel",
    "DAUResult",
    "DagPeptides",
    "GroupingScheme",
    "Proteome",
    "build_background_model",
    "build_background_model_from_samples",
    "build_windows_from_flanks",
    "build_windows_from_positions",
    "format_sequence",
    "load_background_samples",
    "make_daglogo",
    "plot_logo",
    "prepare_proteome",
    "run_daglogo",
    "test_dau",
    "to_logo_matrix",
]
