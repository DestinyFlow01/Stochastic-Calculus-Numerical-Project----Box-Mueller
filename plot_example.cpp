#include <Python.h>
#include <iostream>

void initializePython() {
    Py_Initialize();
}

void finalizePython() {
    Py_Finalize();
}

void createPlot() {
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    PyRun_SimpleString("plt.plot([1, 2, 3, 4], [1, 4, 9, 16])");
    PyRun_SimpleString("plt.ylabel('some numbers')");
    PyRun_SimpleString("plt.show()");
}

int main() {
    initializePython();
    createPlot();
    finalizePython();

    return 0;
}
