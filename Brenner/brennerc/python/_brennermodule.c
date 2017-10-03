/* Python interface to Brenner molecular dynamics library
 *
 * Written by Peter McCluskey
 */

#ifndef Float
#define Float float
#define Double double
#endif
#include "Python.h"
#include "brenner.h"

#ifndef True
#define True 1
#endif

#ifndef False
#define False 0
#endif

typedef struct
{
  const char *dx3dt3;
  const char *accel;
  const char *velocity;
  const char *movable;
  const char *thermostated;
  const char *atom_type;
  int type_is_symbol;
  double cvt_dist_unit;
  const char *starttime;
  const char *timestep;
} FieldMap;

typedef struct {
  PyObject_HEAD
  BrennerMainInfo *info;
  int init_done;
  int kflag;
  const FieldMap *fields;
} PyBrennerObject;

static void
brenner_dealloc(PyBrennerObject *self)
{
  free(self->info);
  PyMem_DEL(self);
}

static PyObject*
brenner_call(PyBrennerObject *self, PyObject *args, PyObject *keywords)
{
  PyObject *obj;
  PyObject *atom_list;
  int steps = 0;
  int num_atoms;
  int i, index;
  const FieldMap *fields = self->fields;
  BrennerMainInfo *info = self->info;
  if (!PyArg_ParseTuple(args, "O|i", &obj, &steps))
    return NULL;
  atom_list = PyObject_CallMethod(obj, "atomList", NULL);
  if(!atom_list)
    return NULL;
  num_atoms = PyList_Size(atom_list);
  if(num_atoms != info->num_atms)
  {
    char buf[128];
    sprintf(buf, "num_atoms in object passed (%d) is not equal to number in original object %d\n",
	    num_atoms, info->num_atms);
    Py_DECREF(atom_list);
    PyErr_SetString(PyExc_ValueError, buf);
    return NULL;
  }
  if(!self->init_done)
  {
    init_bren(self->info, self->kflag, 0);
    self->init_done = 1;
  }
  if(self->kflag != 6)
  {
	int j;
	for (j = 0; j < steps; ++j) {
		bren_1_step(self->info, self->kflag);
	}
  }
  for(index = 0; index < info->num_atms; ++index)
  {
    PyObject *a_py = PyList_GetItem(atom_list, index);
    PyObject *position = PyList_New(3);
    const BrenAtom *aptr = &info->atm_num[index];
    PyObject *x_py = PyFloat_FromDouble(aptr->coord.x / fields->cvt_dist_unit);
    PyObject *y_py = PyFloat_FromDouble(aptr->coord.y / fields->cvt_dist_unit);
    PyObject *z_py = PyFloat_FromDouble(aptr->coord.z / fields->cvt_dist_unit);
    PyList_SetItem(position, 0, x_py);
    PyList_SetItem(position, 1, y_py);
    PyList_SetItem(position, 2, z_py);
    if(!PyObject_CallMethod(a_py, "setPosition", "O", position))
      return NULL;
    Py_DECREF(position);
  }
  Py_DECREF(atom_list);
  Py_INCREF(Py_None);
  return Py_None;
}

static char PyBrenner_Type__doc__[] = 
  "interface to Brenner molecular dynamics C code";


/* Type object */

statichere PyTypeObject PyBrenner_Type = {
  PyObject_HEAD_INIT(NULL)
  0,			          /*ob_size*/
  "Brenner",		          /*tp_name*/
  sizeof(PyBrennerObject),	  /*tp_basicsize*/
  0,			          /*tp_itemsize*/
  /* methods */
  (destructor)brenner_dealloc,   /*tp_dealloc*/
  0,			          /*tp_print*/
  0,/*(getattrfunc)brenner_getattr, */ /*tp_getattr*/
  0, 			          /*tp_setattr*/
  0,			          /*tp_compare*/
  0,                              /*tp_repr*/
  0,                              /*tp_as_number*/
  0,			          /*tp_as_sequence*/
  0,			          /*tp_as_mapping*/
  0,			          /*tp_hash*/
  (ternaryfunc)brenner_call,	  /*tp_call*/
  0,                              /*tp_str*/
  0,                              /*tp_getattro*/
  0,                              /*tp_setattro*/
  /* Space for future expansion */
  0L,0L,
  /* Documentation string */
  PyBrenner_Type__doc__
};


static const FieldMap mmtk_fields = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  "name",
  True,
  10.0,
  "starttime",
  "timestep"
};

static const FieldMap std_fields = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  "atomic_number",
  False,
  1.0,
  "starttime",
  "timestep"
};

/* GetElemNumber is adapted from RasMol 2.6 */

static int GetElemNumber(const char *sym)
{
    char ch1 = toupper(sym[0]), ch2 = toupper(sym[1]);

    switch( ch1 )
    {   case(' '):  switch( ch2 )
                    {   case('B'):  return(  5 );
                        case('C'):  return(  6 );
                        case('D'):  return(  1 );
                        case('F'):  return(  9 );
                        case('H'):  return(  1 );
                        case('I'):  return( 53 );
                        case('K'):  return( 19 );
                        case('L'):  return(  1 );
                        case('N'):  return(  7 );
                        case('O'):  return(  8 );
                        case('P'):  return( 15 );
                        case('S'):  return( 16 );
                        case('U'):  return( 92 );
                        case('V'):  return( 23 );
                        case('W'):  return( 74 );
                        case('Y'):  return( 39 );
                    }
                    break;

        case('A'):  switch( ch2 )
                    {   case('C'):  return( 89 );
                        case('G'):  return( 47 );
                        case('L'):  return( 13 );
                        case('M'):  return( 95 );
                        case('R'):  return( 18 );
                        case('S'):  return( 33 );
                        case('T'):  return( 85 );
                        case('U'):  return( 79 );
                    }
                    break;

        case('B'):  switch( ch2 )
                    {   case('A'):  return( 56 );
                        case('E'):  return(  4 );
                        case('I'):  return( 83 );
                        case('K'):  return( 97 );
                        case('R'):  return( 35 );
                    }
                    break;

        case('C'):  switch( ch2 )
                    {   case('A'):  return( 20 );
                        case('D'):  return( 48 );
                        case('E'):  return( 58 );
                        case('F'):  return( 98 );
                        case('L'):  return( 17 );
                        case('M'):  return( 96 );
                        case('O'):  return( 27 );
                        case('R'):  return( 24 );
                        case('S'):  return( 55 );
                        case('U'):  return( 29 );
                    }
                    break;

        case('D'):  if( ch2=='Y' )
                        return( 66 );
                    break;

        case('E'):  if( ch2=='R' )
                    {   return( 68 );
                    } else if( ch2=='S' )
                    {   return( 99 );
                    } else if( ch2=='U' )
                        return( 63 );
                    break;

        case('F'):  if( ch2=='E' )
                    {   return(  26 );
                    } else if( ch2=='M' )
                    {   return( 100 );
                    } else if( ch2=='R' )
                        return(  87 );
                    break;

        case('G'):  if( ch2=='A' )
                    {   return( 31 );
                    } else if( ch2=='D' )
                    {   return( 64 );
                    } else if( ch2=='E' )
                        return( 32 );
                    break;

        case('H'):  if( ch2=='E' )
                    {   return(  2 );
                    } else if( ch2=='F' )
                    {   return( 72 );
                    } else if( ch2=='G' )
                    {   return( 80 );
                    } else if( ch2=='O' )
                        return( 67 );
                    break;

        case('I'):  if( ch2=='N' )
                    {   return( 49 );
                    } else if( ch2=='R' )
                        return( 77 );
                    break;

        case('K'):  if( ch2=='R' )
                        return( 36 );
                    break;

        case('L'):  if( ch2=='A' )
                    {   return(  57 );
                    } else if( ch2=='I' )
                    {   return(   3 );
                    } else if( (ch2=='R') || (ch2=='W') )
                    {   return( 103 );
                    } else if( ch2=='U' )
                        return(  71 );
                    break;

        case('M'):  if( ch2=='D' )
                    {   return( 101 );
                    } else if( ch2=='G' )
                    {   return(  12 );
                    } else if( ch2=='N' )
                    {   return(  25 );
                    } else if( ch2=='O' )
                        return(  42 );
                    break;

        case('N'):  switch( ch2 )
                    {   case('A'):  return(  11 );
                        case('B'):  return(  41 );
                        case('D'):  return(  60 );
                        case('E'):  return(  10 );
                        case('I'):  return(  28 );
                        case('O'):  return( 102 );
                        case('P'):  return(  93 );
                    }
                    break;

        case('O'):  if( ch2=='S' )
                        return( 76 );
                    break;

        case('P'):  switch( ch2 )
                    {   case('A'):  return( 91 );
                        case('B'):  return( 82 );
                        case('D'):  return( 46 );
                        case('M'):  return( 61 );
                        case('O'):  return( 84 );
                        case('R'):  return( 59 );
                        case('T'):  return( 78 );
                        case('U'):  return( 94 );
                    }
                    break;

        case('R'):  switch( ch2 )
                    {   case('A'):  return( 88 );
                        case('B'):  return( 37 );
                        case('E'):  return( 75 );
                        case('H'):  return( 45 );
                        case('N'):  return( 86 );
                        case('U'):  return( 44 );
                    }
                    break;

        case('S'):  switch( ch2 )
                    {   case('B'):  return( 51 );
                        case('C'):  return( 21 );
                        case('E'):  return( 34 );
                        case('I'):  return( 14 );
                        case('M'):  return( 62 );
                        case('N'):  return( 50 );
                        case('R'):  return( 38 );
                    }
                    break;

        case('T'):  switch( ch2 )
                    {   case('A'):  return( 73 );
                        case('B'):  return( 65 );
                        case('C'):  return( 43 );
                        case('E'):  return( 52 );
                        case('H'):  return( 90 );
                        case('I'):  return( 22 );
                        case('L'):  return( 81 );
                        case('M'):  return( 69 );
                    }
                    break;

        case('X'):  if( ch2=='E' )
                        return( 54 );
                    break;

        case('Y'):  if( ch2=='B' )
                        return( 70 );
                    break;

        case('Z'):  if( ch2=='N' )
                    {   return( 30 );
                    } else if( ch2=='R' )
                        return( 40 );
                    break;
    }

    if( (ch1>='0') && (ch1<='9') )
        if( (ch2=='H') || (ch2=='D') )
            return( 1 ); /* Hydrogen */

    return( 0 );
}

static int
init_coord(BrennerMainInfo *info, PyObject *object, const FieldMap *fields)
{
  BrenAtom *atm_num = info->atm_num;
  int index, natom;
  int num_atoms;
  int dummy;
  double minx = 1.e99, maxx = -1.e99;
  double miny = 1.e99, maxy = -1.e99;
  double minz = 1.e99, maxz = -1.e99;
  PyObject *atom_list;
  PyObject *x, *y, *z;
  PyObject *position;
  PyObject* number0 = PyInt_FromLong(0);
  PyObject* number1 = PyInt_FromLong(1);
  PyObject* number2 = PyInt_FromLong(2);

  if(fields->timestep)
  {
    PyObject *py_ts = PyObject_GetAttrString(object, fields->timestep);
    if(py_ts)
    {
      info->timestep = PyFloat_AsDouble(py_ts);
    }
    else PyErr_Clear();
  }
  if(fields->starttime)
  {
    PyObject *py_st = PyObject_GetAttrString(object, fields->starttime);
    if(py_st)
    {
      info->starttime = PyFloat_AsDouble(py_st);
    }
    else PyErr_Clear();
  }
  atom_list = PyObject_CallMethod(object, "atomList", NULL);
  if(!atom_list)
    return False;
  if(!PyObject_IsTrue(atom_list))
  {
    Py_DECREF(atom_list);
    return False;
  }
  num_atoms = PyList_Size(atom_list);
  if(num_atoms > MAX_ATOMS)
  {
    char buf[128];
    sprintf(buf, "num_atoms %d greater than MAX_ATOMS %d\n",
	    num_atoms, MAX_ATOMS);
    PyErr_SetString(PyExc_ValueError, buf);
    return False;
  }
  for(index = 0; index < num_atoms; ++index)
  {
    PyObject* py_atom_type;
    PyObject* a_py = PyList_GetItem(atom_list, index);
    position = PyObject_CallMethod(a_py, "position", NULL);
    if(!position)
      return False;
    if(!PyObject_IsTrue(position))
    {
      char buf[512];
      sprintf(buf,"Atom #%d has no position; ignored\n",index);
      fprintf(stderr,buf);
      Py_DECREF(position);
      return True;
    }
    x = PyObject_GetItem(position, number0);
    if(!x) return False;
    atm_num[index].coord.x = PyFloat_AsDouble(x) * fields->cvt_dist_unit;
    if(atm_num[index].coord.x > maxx) maxx = atm_num[index].coord.x;
    if(atm_num[index].coord.x < minx) minx = atm_num[index].coord.x;
    y = PyObject_GetItem(position, number1);
    if(!y) return False;
    atm_num[index].coord.y = PyFloat_AsDouble(y) * fields->cvt_dist_unit;
    if(atm_num[index].coord.y > maxy) maxy = atm_num[index].coord.y;
    if(atm_num[index].coord.y < miny) miny = atm_num[index].coord.y;
    z = PyObject_GetItem(position, number2);
    if(!z) return False;
    atm_num[index].coord.z = PyFloat_AsDouble(z) * fields->cvt_dist_unit;
    if(atm_num[index].coord.z > maxz) maxz = atm_num[index].coord.z;
    if(atm_num[index].coord.z < minz) minz = atm_num[index].coord.z;
    Py_DECREF(position);

    py_atom_type = PyObject_GetAttrString(a_py, fields->atom_type);
    if(!py_atom_type) return False;
    if(fields->type_is_symbol)
    {
      const char *sym = PyString_AsString(py_atom_type);
      if(!sym || !*sym)
      {
	char buf[128];
	sprintf(buf, "Atom %d has an empty symbol", index);
	PyErr_SetString(PyExc_ValueError, buf);
	return False;
      }
      natom = GetElemNumber(sym);
      if(!natom)
      {
	char buf[128];
	sprintf(buf, "Atom %d has an invalid symbol: %s", index, sym);
	PyErr_SetString(PyExc_ValueError, buf);
	return False;
      }
    }
    else
    {
      natom = PyInt_AsLong(py_atom_type);
      if(natom < 1)
      {
	char buf[128];
	sprintf(buf, "Atom %d has an invalid element number: %d", index, natom);
	PyErr_SetString(PyExc_ValueError, buf);
	return False;
      }
    }
    Py_DECREF(py_atom_type);

    if(natom > NTYPES || !info->kt[natom])
    {
	char buf[128];
	sprintf(buf, "Atom %d has an unsupported element number: %d", index, natom);
	PyErr_SetString(PyExc_ValueError, buf);
	return False;
    }
    atm_num[index].ktype = info->kt[natom];
    ++info->noa[info->kt[natom]];

    atm_num[index].number = index;
    atm_num[index].type = natom;
    atm_num[index].mass = info->xmass[info->kt[natom]];
    atm_num[index].thermostated = 0;
    if(fields->thermostated)
    {
      PyObject *py_thermostated = PyObject_GetAttrString(a_py, fields->thermostated);
      if(py_thermostated)
	atm_num[index].thermostated = PyInt_AsLong(py_thermostated);
      Py_DECREF(py_thermostated);
    }
    atm_num[index].movable = 1;
    if(fields->movable)
    {
      PyObject *py_movable = PyObject_GetAttrString(a_py, fields->movable);
      if(py_movable)
	atm_num[index].movable = PyInt_AsLong(py_movable);
      Py_DECREF(py_movable);
    }
    if(fields->velocity)
    {
      PyObject *py_vel = PyObject_GetAttrString(a_py, fields->velocity);
      if(!py_vel)
	return False;
      x = PyObject_GetItem(py_vel, number0);
      if(!x) return False;
      atm_num[index].velocity.x = PyFloat_AsDouble(x) * info->timestep;
      y = PyObject_GetItem(py_vel, number1);
      if(!y) return False;
      atm_num[index].velocity.y = PyFloat_AsDouble(y) * info->timestep;
      z = PyObject_GetItem(py_vel, number2);
      if(!z) return False;
      atm_num[index].velocity.z = PyFloat_AsDouble(z) * info->timestep;
      Py_DECREF(py_vel);
    }
#if 0
    else if(atm_num[index].movable)
    {
      const double JITTER = 1.e-3;
      atm_num[index].velocity.x = (rand()/(double)RAND_MAX)*2*JITTER-JITTER;
      atm_num[index].velocity.y = (rand()/(double)RAND_MAX)*2*JITTER-JITTER;
      atm_num[index].velocity.z = (rand()/(double)RAND_MAX)*2*JITTER-JITTER;
    }
#endif
    if(fields->accel)
    {
      PyObject *py_acc = PyObject_GetAttrString(a_py, fields->accel);
      if(!py_acc)
	return False;
      x = PyObject_GetItem(py_acc, number0);
      if(!x) return False;
      atm_num[index].accel.x = PyFloat_AsDouble(x);
      y = PyObject_GetItem(py_acc, number1);
      if(!y) return False;
      atm_num[index].accel.y = PyFloat_AsDouble(y);
      z = PyObject_GetItem(py_acc, number2);
      if(!z) return False;
      atm_num[index].accel.z = PyFloat_AsDouble(z);
      Py_DECREF(py_acc);
    }
    else if(atm_num[index].movable)
    {
      const double small = 1.e-12;
      atm_num[index].accel.x = (rand()/(double)RAND_MAX)*2*small-small;
      atm_num[index].accel.y = (rand()/(double)RAND_MAX)*2*small-small;
      atm_num[index].accel.z = (rand()/(double)RAND_MAX)*2*small-small;
    }
    if(fields->dx3dt3)
    {
      PyObject *py_dx3 = PyObject_GetAttrString(a_py, fields->dx3dt3);
      if(!py_dx3)
	return False;
      x = PyObject_GetItem(py_dx3, number0);
      if(!x) return False;
      atm_num[index].dx3dt3.x = PyFloat_AsDouble(x);
      y = PyObject_GetItem(py_dx3, number1);
      if(!y) return False;
      atm_num[index].dx3dt3.y = PyFloat_AsDouble(y);
      z = PyObject_GetItem(py_dx3, number2);
      if(!z) return False;
      atm_num[index].dx3dt3.z = PyFloat_AsDouble(z);
      Py_DECREF(py_dx3);
    }
    else if(atm_num[index].movable)
    {
      const double small = 1.e-12;
      atm_num[index].dx3dt3.x = (rand()/(double)RAND_MAX)*2*small-small;
      atm_num[index].dx3dt3.y = (rand()/(double)RAND_MAX)*2*small-small;
      atm_num[index].dx3dt3.z = (rand()/(double)RAND_MAX)*2*small-small;
    }
  }
  info->cube[0] = 2*ceil(maxx - minx);
  info->cube[1] = 2*ceil(maxy - miny);
  info->cube[2] = 2*ceil(maxz - minz);
  if(out_of_box(info->cube[0]) || out_of_box(info->cube[1])
     || out_of_box(info->cube[2]))
    printf("warning: large cube size (%.2f,%.2f,%.2f) will require slower algorithm\n",
	   info->cube[0], info->cube[1], info->cube[2]);
  info->num_atms = num_atoms;
  for(index = 0; index < num_atoms; ++index)
  {
    if(atm_num[index].coord.x > info->cube[0]/2)
      atm_num[index].coord.x -= info->cube[0];
    if(atm_num[index].coord.y > info->cube[1]/2)
      atm_num[index].coord.y -= info->cube[1];
    if(atm_num[index].coord.z > info->cube[2]/2)
      atm_num[index].coord.z -= info->cube[2];
    atm_num[index].prev_coord.x = atm_num[index].coord.x;
    atm_num[index].prev_coord.y = atm_num[index].coord.y;
    atm_num[index].prev_coord.z = atm_num[index].coord.z;
  }
  printf("read %d atoms\n", info->num_atms);
  Py_DECREF(atom_list);
  Py_DECREF(number0);
  Py_DECREF(number1);
  Py_DECREF(number2);
}

static PyObject *
construct(PyObject *dummy, PyObject *args, int kflag, const FieldMap *fields)
{
  PyBrennerObject *self;
  PyObject *object;
  if (!PyArg_ParseTuple(args, "O", &object))
    return NULL;
  self = PyObject_NEW(PyBrennerObject, &PyBrenner_Type);
  if (self == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  self->init_done = 0;
  self->info = alloc_bren(&kflag);
  self->kflag = kflag;
  self->fields = fields;

  if(!init_coord(self->info, object, fields))
    return False;
  return (PyObject *)self;
}

static PyObject *
construct_Minimizer(PyObject *dummy, PyObject *args)
{
  return construct(dummy, args, 6, &std_fields);
}

static PyObject *
construct_Integrator(PyObject *dummy, PyObject *args)
{
  return construct(dummy, args, 1, &std_fields);
}

static PyObject *
construct_MMTK_Minimizer(PyObject *dummy, PyObject *args)
{
  return construct(dummy, args, 6, &mmtk_fields);
}

static PyObject *
construct_MMTK_Integrator(PyObject *dummy, PyObject *args)
{
  return construct(dummy, args, 1, &mmtk_fields);
}

/*
 * List of functions defined in the module
 */

static PyMethodDef brenner_methods[] = {
  {"BrennerMinimizer", construct_Minimizer, METH_VARARGS},
  {"BrennerIntegrator", construct_Integrator, METH_VARARGS},
  {"MMTKBrennerMinimizer", construct_MMTK_Minimizer, METH_VARARGS},
  {"MMTKBrennerIntegrator", construct_MMTK_Integrator, METH_VARARGS},
  {NULL, NULL}		/* sentinel */
};

/* Initialization function for the module */

void
#ifdef IBMPC
__declspec(dllexport)
#endif
init_brenner()
{
  PyObject *m;

  /* Create the module and add the functions */
  m = Py_InitModule("_brenner", brenner_methods);
  
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module _brenner");
}

void
my_exit(int code)
{
  exit(code);
}
