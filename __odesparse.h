/* This file should be included in the multipack module */
/* $Revision: 5942 $ */
/*  module_methods:
 {"odeint", (PyCFunction) odesparse_odeints, METH_VARARGS|METH_KEYWORDS, doc_odeints},
 */
/* link libraries: (should be listed in separate lines)
   odesparse
 */
/* python files: (to be appended to Multipack.py)
   odesparse.py
 */
#include <dlfcn.h>
#include <stdio.h>

#if defined(UPPERCASE_FORTRAN)
	#if defined(NO_APPEND_FORTRAN)
		/* nothing to do here */
	#else
		#define DLSODES  DLSODES_
	#endif
#else
	#if defined(NO_APPEND_FORTRAN)
		#define DLSODES  dlsodes
	#else
		#define DLSODES  dlsodes_
	#endif
#endif

void DLSODES();

void ode_function(int *n, double *t, double *y, double *ydot)
{
 /* This is the function called from the Fortran code it should
        -- use call_python_function to get a multiarrayobject result
	-- check for errors and return -1 if any
 	-- otherwise place result of calculation in ydot
  */

  PyArrayObject *result_array = NULL;
  PyObject *arg1, *arglist;

  /* Append t to argument list */
  if ((arg1 = PyTuple_New(1)) == NULL) {
    if (PyErr_Occurred()) {
      PyErr_Print();
    }
    return;
  }
  PyTuple_SET_ITEM(arg1, 0, PyFloat_FromDouble(*t)); 
                /* arg1 now owns newly created reference */
  if ((arglist = PySequence_Concat( arg1, multipack_extra_arguments)) == NULL) {
    if (PyErr_Occurred()) {
      PyErr_Print();
    }
    Py_DECREF(arg1);
    return;
  }
  Py_DECREF(arg1);    /* arglist has reference */
  
  result_array = (PyArrayObject *)call_python_function(multipack_python_function, *n, y, arglist, 1, odesparse_error);
  if (result_array == NULL) {
    PyErr_Print();
    Py_DECREF(arglist);
    return;
  }
  memcpy(ydot, result_array->data, (*n)*sizeof(double));
  Py_DECREF(result_array);
  Py_DECREF(arglist);
  return;
}

int setup_extra_inputs(PyArrayObject **ap_rtol, PyObject *o_rtol, PyArrayObject **ap_atol, PyObject *o_atol, PyArrayObject **ap_tcrit, PyObject *o_tcrit, int *numcrit, int neq)
{
  int itol = 0;
  double tol=1.49012e-8;
  npy_intp one = 1;

  /* Setup tolerances */
  if (o_rtol == NULL) {
    *ap_rtol = (PyArrayObject *)PyArray_SimpleNew(1, &one, PyArray_DOUBLE);
    if (*ap_rtol == NULL) PYERR2(odesparse_error,"Error constructing relative tolerance.");
    *(double *)(*ap_rtol)->data = tol;                /* Default */
  }
  else {
    *ap_rtol = (PyArrayObject *)PyArray_ContiguousFromObject(o_rtol,PyArray_DOUBLE,0,1);
    if (*ap_rtol == NULL) PYERR2(odesparse_error,"Error converting relative tolerance.");
    if ((*ap_rtol)->nd == 0); /* rtol is scalar */
    else if ((*ap_rtol)->dimensions[0] == neq)
      itol |= 2;      /* Set rtol array flag */
    else
      PYERR(odesparse_error,"Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
  }

  if (o_atol == NULL) {
    *ap_atol = (PyArrayObject *)PyArray_SimpleNew(1,&one,PyArray_DOUBLE);
    if (*ap_atol == NULL) PYERR2(odesparse_error,"Error constructing absolute tolerance");
    *(double *)(*ap_atol)->data = tol;
  }
  else {
    *ap_atol = (PyArrayObject *)PyArray_ContiguousFromObject(o_atol,PyArray_DOUBLE,0,1);
    if (*ap_atol == NULL) PYERR2(odesparse_error,"Error converting absolute tolerance.");
    if ((*ap_atol)->nd == 0); /* atol is scalar */
    else if ((*ap_atol)->dimensions[0] == neq) 
      itol |= 1;        /* Set atol array flag */
    else
      PYERR(odesparse_error,"Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
  }
  itol++;             /* increment to get correct value */


  /* Setup t-critical */
  if (o_tcrit != NULL) {
    *ap_tcrit = (PyArrayObject *)PyArray_ContiguousFromObject(o_tcrit,PyArray_DOUBLE,0,1);
    if (*ap_tcrit == NULL) PYERR2(odesparse_error,"Error constructing critical times.");
    *numcrit = PyArray_Size((PyObject *)(*ap_tcrit));
  }
  return itol;

 fail:       /* Needed for use of PYERR */
  return -1;
}


int compute_lrw_liw(int *lrw, int *liw, int neq, int mxords, int nnz, PyObject* o_ia, PyObject* o_ja)
{
  int lwm;
  if (mxords < 0) PYERR(odesparse_error,"Incorrect value for mxords");
  if(o_ia!=NULL && o_ja!=NULL)
    {
      lwm = 2*nnz + 2*neq + (nnz+10*neq);
      *lrw = NPY_MAX(*lrw,20 +  (mxords+4)*neq + lwm);
      *liw = 31 + neq + nnz +100;
    }
  else
    { 
      lwm = 2*nnz + 2*neq + (nnz+10*neq);
      *lrw = NPY_MAX(*lrw,20 +  (mxords+4)*neq + lwm);
      *liw = 30;
    }
  return 0;

 fail:
  return -1;
}


static char doc_odeints[] = "[y,{infodict,}istate] = odeints(fun, y0, t, args=(), full_output=0, rtol=, atol=, tcrit=, h0=0.0, hmax=0.0, hmin=0.0, mxstep=0.0, mxhnil=0, mxords=0)\n  yprime = fun(y,t,...)";

static PyObject *odesparse_odeints(PyObject *dummy, PyObject *args, PyObject *kwdict) {
  PyObject *fcn, *y0, *p_tout, *o_rtol=NULL, *o_atol=NULL, *o_ia=NULL, *o_ja=NULL;
  PyArrayObject *ap_y = NULL, *ap_yout= NULL;
  PyArrayObject *ap_rtol=NULL, *ap_atol=NULL;
  PyArrayObject *ap_tout = NULL;
  PyObject *extra_args = NULL;
  int      neq, itol=1, itask=1, istate=1, iopt=0, lrw=0, *iwork, liw;
  double   *y, t, *tout, *rtol, *atol, *rwork;
  double   h0=0.0, hmax=0.0, hmin=0.0, seth=0.0;
  int      mxstep=0, mxhnil=0, mxords=5, nnz=0;
  PyObject *o_tcrit=NULL;
  PyArrayObject *ap_tcrit=NULL;
  PyArrayObject *ap_hu=NULL, *ap_tcur=NULL, *ap_tolsf=NULL;
  PyArrayObject *ap_ia=NULL, *ap_ja=NULL;
  PyArrayObject *ap_nst=NULL, *ap_nfe=NULL, *ap_nje=NULL, *ap_nqu=NULL;
  int      imxer=0, lenrw=0, leniw=0, ngp = 0, nzl = 0, nzu = 0, mf=222, nlu=0;
  npy_intp out_sz=0,dims[2];
  int      k, ntimes, crit_ind=0;
  int      allocated = 0, full_output = 0, numcrit=0;
  double   *yout, *yout_ptr, *tout_ptr, *tcrit;
  double   *wa;
  static char *kwlist[] = {"fun","y0","t","args","full_output","rtol","atol","tcrit","h0","hmax","hmin","seth", "mxords", "mxstep","mxhnil","lrw", "nnz", "ia", "ja", NULL};
  int i;
  STORE_VARS();


  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "OOO|OiOOOddddiiiiiOO", kwlist, &fcn, &y0, &p_tout, &extra_args, &full_output, &o_rtol, &o_atol, &o_tcrit, &h0, &hmax, &hmin, &seth, &mxords, &mxstep, &mxhnil, &lrw, &nnz, &o_ia, &o_ja)) return NULL;

  if (o_tcrit == Py_None) {
    o_tcrit = NULL;
  }
  if (o_rtol == Py_None) {
    o_rtol = NULL;
  }
  if (o_atol == Py_None) {
    o_atol = NULL;
  }
  if (o_ia == Py_None) {
    o_ia = NULL;
  }
  if (o_ja == Py_None) {
    o_ja = NULL;
  }

  INIT_FUNC(fcn,extra_args,odesparse_error);

  /* Initial input vector */
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y0, PyArray_DOUBLE, 0, 1);
  if (ap_y == NULL) goto fail;
  y = (double *) ap_y->data;
  neq = PyArray_Size((PyObject *)ap_y);
  dims[1] = neq;

  /* Set of output times for integration */
  ap_tout = (PyArrayObject *)PyArray_ContiguousFromObject(p_tout, PyArray_DOUBLE, 0, 1);
  if (ap_tout == NULL) goto fail;
  tout = (double *)ap_tout->data;
  ntimes = PyArray_Size((PyObject *)ap_tout);
  dims[0] = ntimes;
  t = tout[0];

  /* Setup array to hold the output evaluations*/
  ap_yout= (PyArrayObject *)PyArray_SimpleNew(2,dims,PyArray_DOUBLE);
  if (ap_yout== NULL) goto fail;
  yout = (double *) ap_yout->data;
  /* Copy initial vector into first row of output */
  memcpy(yout, y, neq*sizeof(double));  /* copy intial value to output */
  yout_ptr = yout + neq;    /* set output pointer to next position */

  itol = setup_extra_inputs(&ap_rtol, o_rtol, &ap_atol, o_atol, &ap_tcrit, o_tcrit, &numcrit, neq);
  if (itol < 0 ) goto fail;  /* Something didn't work */
  rtol = (double *) ap_rtol->data;
  atol = (double *) ap_atol->data;
  if (o_tcrit != NULL) tcrit = (double *)(ap_tcrit->data);




  /* Find size of working arrays*/
  if (compute_lrw_liw(&lrw, &liw, neq, mxords, nnz, o_ia, o_ja) < 0) goto fail;

  if ((wa = (double *)malloc(lrw*sizeof(double) + liw*sizeof(int)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }

  allocated = 1;
  rwork = wa;
  iwork = (int *)(wa + lrw);


  if(o_ia !=NULL && o_ja!=NULL) {
    ap_ia=(PyArrayObject *) o_ia;
    ap_ja=(PyArrayObject *) o_ja;
    //ap_ia = (PyArrayObject *)PyArray_ContiguousFromObject(o_ia, PyArray_INT, 0, 1);
    //    if (ap_ia == NULL) PYERR2(odesparse_error,"Error constructing ia");
    //ap_ja = (PyArrayObject *)PyArray_ContiguousFromObject(o_ja, PyArray_INT, 0, 1);
    //    if (ap_ja == NULL) PYERR2(odesparse_error,"Error constructing ja");
      for(i=0;i<neq+1;++i) {
	iwork[30+i]=*((int*)(ap_ia->data+i*ap_ia->strides[0]));
      }
      for(i=0;i<nnz;++i) {
	iwork[31+neq+i]=*((int*)(ap_ja->data+i*ap_ja->strides[0]));
      }
      //      printf("array ends: %d %d\n", iwork[31+neq],iwork[31+neq+nnz-1]); 
      mf=22;
  }

  //  printf("mf: %d\n nnz: %d\n", mf, nnz);
  if (h0 != 0.0 || hmax != 0.0 || hmin != 0.0 || seth!=0.0 || mxstep != 0 || mxhnil != 0 || mxords != 0) {
    rwork[4] = h0; rwork[5] = hmax; rwork[6] = hmin; rwork[7]=seth;
    iwork[4] = mxords; iwork[5] = mxstep; iwork[6] = mxhnil;
    iopt = 1;
  }
  istate = 1;
  k = 1;

  /* If full output make some useful output arrays */
  if (full_output) {
    out_sz = ntimes-1;
    ap_hu = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_DOUBLE);
    ap_tcur = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_DOUBLE);
    ap_tolsf = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_DOUBLE);
    ap_nst = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    ap_nfe = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    ap_nje = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    ap_nqu = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    if (ap_hu == NULL || ap_tcur == NULL || ap_tolsf == NULL || ap_nst == NULL || ap_nfe == NULL || ap_nje == NULL || ap_nqu == NULL) goto fail;
  }

  if (o_tcrit != NULL) {itask = 4; rwork[0] = *tcrit;}  /* There are critical points */
  while (k < ntimes && istate > 0) {    /* loop over desired times */

    tout_ptr = tout + k;
    /* Use tcrit if relevant */
    if (itask == 4 && *tout_ptr > *(tcrit + crit_ind)) {crit_ind++; rwork[0] = *(tcrit+crit_ind);}
    if (crit_ind >= numcrit) itask = 1;  /* No more critical values */

    DLSODES(ode_function, &neq, y, &t, tout_ptr, &itol, rtol, atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, NULL, &mf);
    if (full_output) {
      *((double *)ap_hu->data + (k-1)) = rwork[10];
      *((double *)ap_tcur->data + (k-1)) = rwork[12];
      *((double *)ap_tolsf->data + (k-1)) = rwork[13];
      *((int *)ap_nst->data + (k-1)) = iwork[10];
      *((int *)ap_nfe->data + (k-1)) = iwork[11];
      *((int *)ap_nje->data + (k-1)) = iwork[12];
      *((int *)ap_nqu->data + (k-1)) = iwork[13];

      if (istate == -5 || istate == -4) {
        imxer = iwork[15];
      } else {
        imxer = -1;
      }
      lenrw = iwork[16];
      leniw = iwork[17];
      nnz = iwork[18];
      ngp = iwork[19];
      nlu = iwork[20];
      nzl = iwork[24];
      nzu = iwork[25];
    }
    if (PyErr_Occurred()) goto fail;
    memcpy(yout_ptr, y, neq*sizeof(double));  /* copy integration result to output*/
    yout_ptr += neq;  k++;
  }

  RESTORE_FUNC();
  Py_DECREF(extra_args);
  Py_DECREF(ap_atol);
  Py_DECREF(ap_rtol);
  Py_XDECREF(ap_tcrit);
  Py_DECREF(ap_y);
  Py_DECREF(ap_tout);
  // Py_XDECREF(ap_ia);
  // Py_XDECREF(ap_ja);
  free(wa);

  /* Do Full output */
    if (full_output) {
      return Py_BuildValue("N{s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i}i",PyArray_Return(ap_yout),
                      "hu",PyArray_Return(ap_hu),
                      "tcur",PyArray_Return(ap_tcur),
                      "tolsf",PyArray_Return(ap_tolsf),
                      "nst",PyArray_Return(ap_nst),
                      "nfe",PyArray_Return(ap_nfe),
                      "nje",PyArray_Return(ap_nje),
                      "nqu",PyArray_Return(ap_nqu),
                      "imxer",imxer,
                      "lenrw",lenrw,
                      "leniw",leniw,
		      "nnz", nnz,
		      "ngp", ngp,
		      "nlu", nlu,
		      "nzl", nzl,
		      "nzu", nzu,
                      istate);
    }
    else {
      return Py_BuildValue("Ni",PyArray_Return(ap_yout),istate);
    }
    
 fail:
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_rtol);
  Py_XDECREF(ap_atol);
  Py_XDECREF(ap_tcrit);
  Py_XDECREF(ap_tout);
  Py_XDECREF(ap_yout);
  //  Py_XDECREF(ap_ia);
  //  Py_XDECREF(ap_ja);
  if (allocated) free(wa);
  if (full_output) {
    Py_XDECREF(ap_hu);
    Py_XDECREF(ap_tcur);
    Py_XDECREF(ap_tolsf);
    Py_XDECREF(ap_nst);
    Py_XDECREF(ap_nfe);
    Py_XDECREF(ap_nje);
    Py_XDECREF(ap_nqu);
  }
  return NULL;
  
}



static char doc_odeints_cfunc[] = "[y,{infodict,}istate] = odeints_cfunc(cfun, y0, t, full_output=0, rtol=, atol=, tcrit=, h0=0.0, hmax=0.0, hmin=0.0, mxstep=0.0, mxhnil=0, mxords=0)\n  yprime = fun(y,t,...)";

static PyObject *odesparse_odeints_cfunc(PyObject *dummy, PyObject *args, PyObject *kwdict) {
  PyObject *y0, *p_tout, *o_rtol=NULL, *o_atol=NULL, *o_ia=NULL, *o_ja=NULL, *o_neq = NULL;
  void *fcn;
  PyArrayObject *ap_neq = NULL;
  PyArrayObject *ap_y = NULL, *ap_yout= NULL;
  PyArrayObject *ap_rtol=NULL, *ap_atol=NULL;
  PyArrayObject *ap_tout = NULL;
  int      neq, ny, itol=1, itask=1, istate=1, iopt=0, lrw=0, *iwork, liw;
  int      *neq_new;
  double   *y, t, *tout, *rtol, *atol, *rwork;
  double   h0=0.0, hmax=0.0, hmin=0.0, seth=0.0;
  int      mxstep=0, mxhnil=0, mxords=5, nnz=0;
  PyObject *o_tcrit=NULL;
  PyArrayObject *ap_tcrit=NULL;
  PyArrayObject *ap_hu=NULL, *ap_tcur=NULL, *ap_tolsf=NULL;
  PyArrayObject *ap_ia=NULL, *ap_ja=NULL;
  PyArrayObject *ap_nst=NULL, *ap_nfe=NULL, *ap_nje=NULL, *ap_nqu=NULL;
  int      imxer=0, lenrw=0, leniw=0, ngp = 0, nzl = 0, nzu = 0, mf=222, nlu=0;
  npy_intp out_sz=0,dims[2];
  int      k, ntimes, crit_ind=0;
  int      allocated = 0, full_output = 0, numcrit=0;
  double   *yout, *yout_ptr, *tout_ptr, *tcrit;
  double   *wa;
  static char *kwlist[] = {"fun","y0","t", "neq", "full_output","rtol","atol","tcrit","h0","hmax","hmin","seth", "mxords", "mxstep","mxhnil","lrw", "nnz", "ia", "ja", NULL};
  int i;
  //  STORE_VARS();





  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "kOO|OiOOOddddiiiiiOO", kwlist, &fcn, &y0, &p_tout, &o_neq, &full_output, &o_rtol, &o_atol, &o_tcrit, &h0, &hmax, &hmin, &seth, &mxords, &mxstep, &mxhnil, &lrw, &nnz, &o_ia, &o_ja)) return NULL;

  //  fprintf(stderr,"fun: %lu\n", (unsigned long) fcn);
  
  //  void* handle = dlopen("/tmp/react.so", RTLD_LAZY);
  //  typedef void (*hello_t)();
  //  hello_t hello = (hello_t) dlsym(handle, "deriv");
  //  fprintf(stderr,"hello: %p %p %ld\n", handle, hello, (long int) hello);
  //  fcn=hello;
 

  if (o_tcrit == Py_None) {
    o_tcrit = NULL;
  }
  if (o_rtol == Py_None) {
    o_rtol = NULL;
  }
  if (o_atol == Py_None) {
    o_atol = NULL;
  }
  if (o_ia == Py_None) {
    o_ia = NULL;
  }
  if (o_ja == Py_None) {
    o_ja = NULL;
  }

  //  INIT_FUNC(fcn,extra_args,odesparse_error);

  /* Initial input vector */
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y0, PyArray_DOUBLE, 0, 1);
  if (ap_y == NULL) goto fail;
  y = (double *) ap_y->data;

  if(o_neq == NULL) goto fail;
  ap_neq = (PyArrayObject *)PyArray_ContiguousFromObject(o_neq, PyArray_INT, 0, 1);
  if(ap_neq == NULL) goto fail;

  neq_new = (int *) ap_neq->data;

  neq = neq_new[0];

  ny = PyArray_Size((PyObject *)ap_y);

  dims[1] = neq;

  /* Set of output times for integration */
  ap_tout = (PyArrayObject *)PyArray_ContiguousFromObject(p_tout, PyArray_DOUBLE, 0, 1);
  if (ap_tout == NULL) goto fail;
  tout = (double *)ap_tout->data;
  ntimes = PyArray_Size((PyObject *)ap_tout);
  dims[0] = ntimes;
  t = tout[0];

  /* Setup array to hold the output evaluations*/
  ap_yout= (PyArrayObject *)PyArray_SimpleNew(2,dims,PyArray_DOUBLE);
  if (ap_yout== NULL) goto fail;
  yout = (double *) ap_yout->data;
  /* Copy initial vector into first row of output */
  memcpy(yout, y, neq*sizeof(double));  /* copy intial value to output */
  yout_ptr = yout + neq;    /* set output pointer to next position */

  itol = setup_extra_inputs(&ap_rtol, o_rtol, &ap_atol, o_atol, &ap_tcrit, o_tcrit, &numcrit, neq);
  if (itol < 0 ) goto fail;  /* Something didn't work */
  rtol = (double *) ap_rtol->data;
  atol = (double *) ap_atol->data;
  if (o_tcrit != NULL) tcrit = (double *)(ap_tcrit->data);

  printf("NEQ, NNZ,: %d %d\n", neq, nnz);

  /* Find size of working arrays*/
  if (compute_lrw_liw(&lrw, &liw, neq, mxords, nnz, o_ia, o_ja) < 0) goto fail;

  if ((wa = (double *)malloc(lrw*sizeof(double) + liw*sizeof(int)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }

  allocated = 1;
  rwork = wa;
  iwork = (int *)(wa + lrw);
     
  if(o_ia !=NULL && o_ja!=NULL) {
    ap_ia=(PyArrayObject *) o_ia;
    ap_ja=(PyArrayObject *) o_ja;
    //   ap_ia = (PyArrayObject *)PyArray_ContiguousFromObject(o_ia, PyArray_INT, 0, 1);
    // if (ap_ia == NULL) PYERR2(odesparse_error,"Error constructing ia");
    //ap_ja = (PyArrayObject *)PyArray_ContiguousFromObject(o_ja, PyArray_INT, 0, 1);
    // if (ap_ja == NULL) PYERR2(odesparse_error,"Error constructing ja");
      for(i=0;i<neq+1;++i) {
	iwork[30+i]=*((int*)(ap_ia->data+i*ap_ia->strides[0]));
      }
      for(i=0;i<nnz;++i) {
	iwork[31+neq+i]=*((int*)(ap_ja->data+i*ap_ja->strides[0]));
      }
      //      printf("array ends: %d %d\n", iwork[31+neq],iwork[31+neq+nnz-1]); 
      mf=22;
  }

  printf("mf: %d\n nnz: %d\n", mf, nnz);
  if (h0 != 0.0 || hmax != 0.0 || hmin != 0.0 || seth!=0.0 || mxstep != 0 || mxhnil != 0 || mxords != 0) {
    rwork[4] = h0; rwork[5] = hmax; rwork[6] = hmin; rwork[7]=seth;
    iwork[4] = mxords; iwork[5] = mxstep; iwork[6] = mxhnil;
    iopt = 1;
  }
  istate = 1;
  k = 1;

  /* If full output make some useful output arrays */
  if (full_output) {
    out_sz = ntimes-1;
    ap_hu = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_DOUBLE);
    ap_tcur = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_DOUBLE);
    ap_tolsf = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_DOUBLE);
    ap_nst = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    ap_nfe = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    ap_nje = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    ap_nqu = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,PyArray_INT);
    if (ap_hu == NULL || ap_tcur == NULL || ap_tolsf == NULL || ap_nst == NULL || ap_nfe == NULL || ap_nje == NULL || ap_nqu == NULL) goto fail;
  }

  if (o_tcrit != NULL) {itask = 4; rwork[0] = *tcrit;}  /* There are critical points */
  while (k < ntimes && istate > 0) {    /* loop over desired times */

    tout_ptr = tout + k;
    /* Use tcrit if relevant */
    if (itask == 4 && *tout_ptr > *(tcrit + crit_ind)) {crit_ind++; rwork[0] = *(tcrit+crit_ind);}
    if (crit_ind >= numcrit) itask = 1;  /* No more critical values */


    DLSODES(fcn, neq_new, y, &t, tout_ptr, &itol, rtol, atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, NULL, &mf);
    if (full_output) {
      *((double *)ap_hu->data + (k-1)) = rwork[10];
      *((double *)ap_tcur->data + (k-1)) = rwork[12];
      *((double *)ap_tolsf->data + (k-1)) = rwork[13];
      *((int *)ap_nst->data + (k-1)) = iwork[10];
      *((int *)ap_nfe->data + (k-1)) = iwork[11];
      *((int *)ap_nje->data + (k-1)) = iwork[12];
      *((int *)ap_nqu->data + (k-1)) = iwork[13];


      if (istate == -5 || istate == -4) {
        imxer = iwork[15];
      } else {
        imxer = -1;
      }
      lenrw = iwork[16];
      leniw = iwork[17];
      nnz = iwork[18];
      ngp = iwork[19];
      nlu = iwork[20];
      nzl = iwork[24];
      nzu = iwork[25];
    }
    if (PyErr_Occurred()) goto fail;
    memcpy(yout_ptr, y, neq*sizeof(double));  /* copy integration result to output*/
    yout_ptr += neq;  k++;
  }

  //  RESTORE_FUNC();
  //  Py_DECREF(extra_args);
  Py_XDECREF(ap_neq);
  Py_DECREF(ap_atol);
  Py_DECREF(ap_rtol);
  Py_XDECREF(ap_tcrit);
  Py_DECREF(ap_y);
  Py_DECREF(ap_tout);
  // Py_XDECREF(ap_ia);
  // Py_XDECREF(ap_ja);
  free(wa);

  /* Do Full output */
    if (full_output) {
      return Py_BuildValue("N{s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i}i",PyArray_Return(ap_yout),
                      "hu",PyArray_Return(ap_hu),
                      "tcur",PyArray_Return(ap_tcur),
                      "tolsf",PyArray_Return(ap_tolsf),
                      "nst",PyArray_Return(ap_nst),
                      "nfe",PyArray_Return(ap_nfe),
                      "nje",PyArray_Return(ap_nje),
                      "nqu",PyArray_Return(ap_nqu),
                      "imxer",imxer,
                      "lenrw",lenrw,
                      "leniw",leniw,
		      "nnz", nnz,
		      "ngp", ngp,
		      "nlu", nlu,
		      "nzl", nzl,
		      "nzu", nzu,
                      istate);
    }
    else {
      return Py_BuildValue("Ni",PyArray_Return(ap_yout),istate);
    }
    
 fail:
    //  RESTORE_FUNC();
  //  Py_XDECREF(extra_args);
  Py_XDECREF(ap_neq);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_rtol);
  Py_XDECREF(ap_atol);
  Py_XDECREF(ap_tcrit);
  Py_XDECREF(ap_tout);
  Py_XDECREF(ap_yout);
  //  Py_XDECREF(ap_ia);
  //  Py_XDECREF(ap_ja);
  if (allocated) free(wa);
  if (full_output) {
    Py_XDECREF(ap_hu);
    Py_XDECREF(ap_tcur);
    Py_XDECREF(ap_tolsf);
    Py_XDECREF(ap_nst);
    Py_XDECREF(ap_nfe);
    Py_XDECREF(ap_nje);
    Py_XDECREF(ap_nqu);
  }
  return NULL;
  
}

static char doc_dl[] = "h = dlopen(dllname,fname)";

void* handle;

static PyObject *odesparse_dlopen(PyObject *self, PyObject *args)
{
  
  const char *dllname, *fname;
  if (!PyArg_ParseTuple(args, "ss", &dllname, &fname))
        return NULL;
  // fprintf(stderr, " dll %s fn %s \n", dllname, fname);

  handle = dlopen(dllname, RTLD_LAZY);
    
  if (!handle) {
    fprintf(stderr, "Cannot open library %s\n", dllname);
        return NULL;
    }
  dlerror();
  void* fn = dlsym(handle, fname);
  const char *dlsym_error = dlerror();
  
  if (dlsym_error) {
    fprintf(stderr, "Cannot load symbol 'hello': %s\n", dlsym_error);
        dlclose(handle);
        return NULL;
  }
  //  fprintf(stderr, " handle %p fn %p \n", handle, fn);
  return Py_BuildValue("k", (long) fn);
}

static PyObject *odesparse_dlclose(PyObject *self, PyObject *args)
{
  if(handle) {
    //    fprintf(stderr, "dllclose %p\n", handle);
    dlclose(handle);
  }
  Py_RETURN_NONE;
}
