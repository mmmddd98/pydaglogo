from __future__ import annotations

import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = DOCS_DIR.parent
SRC_DIR = PACKAGE_ROOT / "src"

sys.path.insert(0, str(SRC_DIR))

project = "pydaglogo"
author = "Markus Kramer and collaborators"
release = "1.0.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path: list[str] = []
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = False
napoleon_numpy_docstring = True

html_theme = "alabaster"
html_static_path: list[str] = []
html_title = "pydaglogo documentation"
