"""
fly-stocker-v2: Modular Drosophila Stock Processing Pipeline

A refactored, modular codebase for processing Drosophila stock data from
FlyBase and BDSC, then organizing split outputs and validating Ref++ stocks.

Pipelines:
    - pipeline_references: Gene → BDSC Stocks → References
    - pipeline_split: Organize stocks by JSON rules and run Ref++ post-split validation

Usage:
    from fly_stocker_v2 import AlleleReferencesPipeline, StockSplitterPipeline
    from fly_stocker_v2.config import Settings
    
    # Pipeline 1: Get allele references
    settings = Settings()
    pipeline = AlleleReferencesPipeline(settings)
    pipeline.run(input_dir, keywords=["sleep", "circadian"], run_functional_validation=False)
    
    # Pipeline 2: Split and validate Ref++ output-sheet stocks
    splitter = StockSplitterPipeline(settings)
    splitter.run(input_dir)
"""

__version__ = "2.0.0"

# Lazy imports to avoid circular dependencies
def __getattr__(name):
    if name == "AlleleReferencesPipeline":
        from .pipeline_references import AlleleReferencesPipeline
        return AlleleReferencesPipeline
    elif name == "StockSplitterPipeline":
        from .pipeline_split import StockSplitterPipeline
        return StockSplitterPipeline
    elif name == "Settings":
        from .config import Settings
        return Settings
    elif name == "is_portable_mode":
        from .config import is_portable_mode
        return is_portable_mode
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "AlleleReferencesPipeline",
    "StockSplitterPipeline",
    "Settings",
    "is_portable_mode",
    "__version__",
]
