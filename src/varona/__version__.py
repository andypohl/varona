import pathlib

__version__ = (
    pathlib.Path(__file__).with_name("__VERSION__").read_text(encoding="utf-8").strip()
)
