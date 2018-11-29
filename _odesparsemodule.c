/*
    Multipack project.
 */
#include "multipack.h"
static PyObject *odesparse_error;
#include "__odesparse.h"

static struct PyMethodDef odesparse_module_methods[] = {
  {"odeints", (PyCFunction) odesparse_odeints, 
   METH_VARARGS|METH_KEYWORDS, doc_odeints},
  {"odeints_cfunc", (PyCFunction) odesparse_odeints_cfunc, 
   METH_VARARGS|METH_KEYWORDS, doc_odeints_cfunc},
  {"dlopen", (PyCFunction) odesparse_dlopen, 
   METH_VARARGS, doc_dl},
  {"dlclose", (PyCFunction) odesparse_dlclose, 
   METH_VARARGS, doc_dl},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_odesparse(void) {
  PyObject *m, *d, *s;
  m = Py_InitModule("_odesparse", odesparse_module_methods);
  import_array();
  d = PyModule_GetDict(m);
  s = PyString_FromString(" 1.9 ");
  PyDict_SetItemString(d, "__version__", s);
  odesparse_error = PyErr_NewException ("odesparse.error", NULL, NULL);
  Py_DECREF(s);
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module odesparse");
}
        
