// -*- C++ -*-
//===---------------------- py_solvers_module.cpp --------------------------===//
//
// Thin CPython + NumPy C-API wrapper around csrc/solvers.cpp (the C++ port
// of aocseq's original Rcpp `solvers.cpp`). This plays the same role the
// Rcpp `// [[Rcpp::export]]` markers played in R: exposing
// GetMahalanobis()/IsoForest() to the host language. Only the Python/NumPy
// C API is used (no pybind11/Boost.Python), to keep the build-time
// dependency footprint to "a C++ compiler + numpy" -- numpy is already a
// runtime dependency of paocseq, and its headers ship with the numpy wheel.
//
//===----------------------------------------------------------------------===//

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <algorithm>
#include <random>
#include <vector>

#include "isolation_forest.hpp"
#include "mahalanobis.hpp"

namespace {

// Convert a 2-D array-like Python object into a contiguous C-contiguous
// float32 numpy array, genes x cells (no copy if it's already exactly
// that dtype/layout).
PyArrayObject *AsFloat32Matrix(PyObject *obj, const char *argname) {
  PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(
      PyArray_FROM_OTF(obj, NPY_FLOAT32, NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST));
  if (arr == nullptr) {
    PyErr_Format(PyExc_TypeError, "%s must be convertible to a 2-D float array", argname);
    return nullptr;
  }
  if (PyArray_NDIM(arr) != 2) {
    PyErr_Format(PyExc_ValueError, "%s must be 2-D (genes x cells)", argname);
    Py_DECREF(arr);
    return nullptr;
  }
  return arr;
}

PyObject *PyMahalanobisDistance(PyObject *self, PyObject *args, PyObject *kwargs) {
  static const char *kwlist[] = {"reference", "query", nullptr};
  PyObject *reference_obj = nullptr;
  PyObject *query_obj = nullptr;

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", const_cast<char **>(kwlist), &reference_obj, &query_obj)) {
    return nullptr;
  }

  PyArrayObject *reference_arr = AsFloat32Matrix(reference_obj, "reference");
  if (!reference_arr) return nullptr;
  PyArrayObject *query_arr = AsFloat32Matrix(query_obj, "query");
  if (!query_arr) {
    Py_DECREF(reference_arr);
    return nullptr;
  }

  npy_intp n_genes = PyArray_DIM(reference_arr, 0);
  npy_intp n_ref = PyArray_DIM(reference_arr, 1);
  npy_intp n_genes_q = PyArray_DIM(query_arr, 0);
  npy_intp n_query = PyArray_DIM(query_arr, 1);

  if (n_genes != n_genes_q) {
    PyErr_SetString(PyExc_ValueError, "reference and query must have the same number of genes (axis 0)");
    Py_DECREF(reference_arr);
    Py_DECREF(query_arr);
    return nullptr;
  }

  const float *reference_data = static_cast<const float *>(PyArray_DATA(reference_arr));
  const float *query_data = static_cast<const float *>(PyArray_DATA(query_arr));

  std::vector<float> reference_vec(reference_data, reference_data + n_genes * n_ref);
  std::vector<float> query_vec(query_data, query_data + n_genes * n_query);

  std::vector<float> distances;
  try {
    distances = paocseq::GetMahalanobis(static_cast<int>(n_genes), static_cast<int>(n_ref),
                                         static_cast<int>(n_query), reference_vec, query_vec);
  } catch (const std::exception &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    Py_DECREF(reference_arr);
    Py_DECREF(query_arr);
    return nullptr;
  }

  Py_DECREF(reference_arr);
  Py_DECREF(query_arr);

  npy_intp out_dims[1] = {n_query};
  PyObject *out = PyArray_SimpleNew(1, out_dims, NPY_FLOAT32);
  if (!out) return nullptr;
  std::copy(distances.begin(), distances.end(), static_cast<float *>(PyArray_DATA(reinterpret_cast<PyArrayObject *>(out))));
  return out;
}

PyObject *PyIsolationForestScore(PyObject *self, PyObject *args, PyObject *kwargs) {
  static const char *kwlist[] = {"reference", "query", "n_trees", "max_height", "random_state", nullptr};
  PyObject *reference_obj = nullptr;
  PyObject *query_obj = nullptr;
  int n_trees = 10;
  int max_height = 30;
  long random_state = 0;
  int have_seed = 0;

  PyObject *random_state_obj = Py_None;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiO", const_cast<char **>(kwlist), &reference_obj, &query_obj,
                                    &n_trees, &max_height, &random_state_obj)) {
    return nullptr;
  }
  if (random_state_obj != Py_None) {
    random_state = PyLong_AsLong(random_state_obj);
    if (random_state == -1 && PyErr_Occurred()) return nullptr;
    have_seed = 1;
  }

  PyArrayObject *reference_arr = AsFloat32Matrix(reference_obj, "reference");
  if (!reference_arr) return nullptr;
  PyArrayObject *query_arr = AsFloat32Matrix(query_obj, "query");
  if (!query_arr) {
    Py_DECREF(reference_arr);
    return nullptr;
  }

  npy_intp n_genes = PyArray_DIM(reference_arr, 0);
  npy_intp n_ref = PyArray_DIM(reference_arr, 1);
  npy_intp n_genes_q = PyArray_DIM(query_arr, 0);
  npy_intp n_query = PyArray_DIM(query_arr, 1);

  if (n_genes != n_genes_q) {
    PyErr_SetString(PyExc_ValueError, "reference and query must have the same number of genes (axis 0)");
    Py_DECREF(reference_arr);
    Py_DECREF(query_arr);
    return nullptr;
  }

  const float *reference_data = static_cast<const float *>(PyArray_DATA(reference_arr));
  const float *query_data = static_cast<const float *>(PyArray_DATA(query_arr));

  std::vector<float> reference_vec(reference_data, reference_data + n_genes * n_ref);
  std::vector<float> query_vec(query_data, query_data + n_genes * n_query);

  unsigned int seed = have_seed ? static_cast<unsigned int>(random_state) : std::random_device{}();

  std::vector<float> scores;
  try {
    scores = paocseq::IsoForest(n_trees, static_cast<int>(n_genes), static_cast<int>(n_ref),
                                 static_cast<int>(n_query), reference_vec, query_vec, max_height, seed);
  } catch (const std::exception &e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    Py_DECREF(reference_arr);
    Py_DECREF(query_arr);
    return nullptr;
  }

  Py_DECREF(reference_arr);
  Py_DECREF(query_arr);

  npy_intp out_dims[1] = {n_query};
  PyObject *out = PyArray_SimpleNew(1, out_dims, NPY_FLOAT32);
  if (!out) return nullptr;
  std::copy(scores.begin(), scores.end(), static_cast<float *>(PyArray_DATA(reinterpret_cast<PyArrayObject *>(out))));
  return out;
}

PyMethodDef kMethods[] = {
    {"mahalanobis_distance", reinterpret_cast<PyCFunction>(PyMahalanobisDistance), METH_VARARGS | METH_KEYWORDS,
     "mahalanobis_distance(reference, query) -> ndarray\n\n"
     "Leave-one-out Mahalanobis distance of each query cell (genes x n_query)\n"
     "to a reference panel (genes x n_ref). Direct C++ port of GetMahalanobis\n"
     "from the original aocseq solvers.cpp."},
    {"isolation_forest_score", reinterpret_cast<PyCFunction>(PyIsolationForestScore), METH_VARARGS | METH_KEYWORDS,
     "isolation_forest_score(reference, query, n_trees=10, max_height=30, random_state=None) -> ndarray\n\n"
     "Kurtosis-driven isolation-forest outlier score of each query cell (genes\n"
     "x n_query) relative to a reference panel (genes x n_ref). Direct C++\n"
     "port of tree()/isoForest() from the original aocseq solvers.cpp."},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef kModuleDef = {
    PyModuleDef_HEAD_INIT,
    "_solvers_ext",
    "Compiled C++ solvers for paocseq (Mahalanobis distance, isolation forest), "
    "ported from the original aocseq R package's solvers.cpp.",
    -1,
    kMethods,
};

} // namespace

PyMODINIT_FUNC PyInit__solvers_ext(void) {
  import_array();
  return PyModule_Create(&kModuleDef);
}
