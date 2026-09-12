"""Build script for paocseq's compiled C++ extension (``paocseq._solvers_ext``).

paocseq is otherwise a pure-Python package (see ``pyproject.toml``); this
``setup.py`` only exists to declare the one ``ext_modules`` entry needed to
compile ``csrc/mahalanobis.cpp`` and ``csrc/isolation_forest.cpp`` (direct
C++ ports of the corrected ``GetMahalanobis.cpp`` and the legacy
``solvers.cpp``'s isolation forest, kept as two separate source files, the
same way the R package's Rcpp code is split to avoid build issues) into
``paocseq._solvers_ext``. If no C++ compiler is available at install time,
``paocseq`` still installs and works, falling back to the pure-Python
implementation in ``paocseq._solvers_py`` (see ``paocseq/_solvers.py``).
"""

from __future__ import annotations

import sys

import numpy as np
from setuptools import Extension, setup

if sys.platform == "win32":
    extra_compile_args = ["/std:c++14", "/O2"]
else:
    extra_compile_args = ["-std=c++14", "-O3"]

ext_modules = [
    Extension(
        name="paocseq._solvers_ext",
        sources=[
            "csrc/mahalanobis.cpp",
            "csrc/isolation_forest.cpp",
            "csrc/py_solvers_module.cpp",
        ],
        include_dirs=[np.get_include(), "csrc"],
        language="c++",
        extra_compile_args=extra_compile_args,
        optional=True,  # don't fail the whole install if this can't compile
    )
]

setup(ext_modules=ext_modules)
