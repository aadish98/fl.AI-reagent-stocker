"""
Canonical package for the fly stocker pipelines.
"""

from .config import Settings, is_portable_mode
from .pipelines.stock_finding import StockFindingPipeline
from .pipelines.stock_splitting import StockSplittingPipeline

__version__ = "3.0.0"

__all__ = [
    "Settings",
    "StockFindingPipeline",
    "StockSplittingPipeline",
    "is_portable_mode",
    "__version__",
]
