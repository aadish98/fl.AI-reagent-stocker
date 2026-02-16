"""Compatibility entry point for `python -m fly_stocker_v2`."""

from .cli import main


if __name__ == "__main__":
    raise SystemExit(main())

