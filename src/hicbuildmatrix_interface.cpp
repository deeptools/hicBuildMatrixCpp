#include "Python.h"
#include "hicbuildmatrix.hpp"

static PyObject* buildMatrix(PyObject* self, PyObject* args) {
    char* forwardReadChar;
    char* reverseReadChar;

    if (!PyArg_ParseTuple(args, "ss", &forwardReadChar, &reverseReadChar))
        return NULL;

    std::string forwardReadString = forwardReadChar;
    std::string reverseReadString = reverseReadChar;


    HiCBuildMatrix* hiCBuildMatrix;
    hiCBuildMatrix = new HiCBuildMatrix(forwardReadString, reverseReadString);
    std:: vector <std::string> pRefId2name;
    //pRefId2name.push_back(10); 
    //pRefId2name.push_back(20); 
    int buffer_size = 100000;
    bool dublication_check = false;
    hiCBuildMatrix->readBamFile(buffer_size, dublication_check, pRefId2name);
    std::string pchrom1 = "chr1";
    int pstart1 = 1000;
    std::string pchrom2 = "chr2";
    int pstart2 = 20000;
    hiCBuildMatrix->is_duplicated( pchrom1, pstart1, pchrom2, pstart2);
    // hiCBuildMatrix->createInitialStructures(forwardReadString);
    return Py_BuildValue("");
}

// definition of available functions for python and which function parsing function in c++ should be called.
static PyMethodDef hicBuildMatrixFunctions[] = {
    {"buildMatrix", (PyCFunction) buildMatrix, METH_VARARGS, "Builds a Hi-C interaction matrix."},
   
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef hicBuildMatrixCppModule = {
    PyModuleDef_HEAD_INIT,
    "hicBuildMatrixCpp",
    NULL,
    -1,
    hicBuildMatrixFunctions
};

// definition of the module for python
PyMODINIT_FUNC PyInit_hicBuildMatrixCpp(void)
{
    return PyModule_Create(&hicBuildMatrixCppModule);
}
