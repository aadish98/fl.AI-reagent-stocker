"""
Entry point for running fly-stocker-v2 as a module.

Usage:
    python -m fly_stocker_v2 get-allele-refs /path/to/input --keywords sleep
    python -m fly_stocker_v2 split-stocks /path/to/input --stocks-per-gene 5
"""

import sys
from .cli import main

if __name__ == "__main__":
    sys.exit(main())
