#include <Python.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace std;

vector<string> split(const string& line, char delimiter) {
	vector<string> result;
	stringstream ss(line);
	string item;
	
	while(getline(ss, item, delimiter)) {
		result.push_back(item);
	}
	return result;
}

void ReadCSV(string filename, vector<double>& Z1, vector<double>& Z2) {
	ifstream file(filename);
    string line;
    bool first_line = true;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            // Skip the header line
            if (first_line) {
                first_line = false;
                continue;
            }

            std::vector<std::string> row = split(line, ',');

            // Convert Z1 and Z2 values from string to double and store them in vectors
            Z1.push_back(std::stod(row[1]));
            Z2.push_back(std::stod(row[2]));
        }
        file.close();
    } 
	else {
        std::cerr << "Unable to open file" << std::endl;
	}
}

void call_Python(vector<double> Z1, vector<double> Z2, vector<double> U1, vector<double> U2) {
    // Convert C++ vectors to Python lists
    PyObject* py_Z1 = PyList_New(Z1.size());
    PyObject* py_Z2 = PyList_New(Z2.size());
    PyObject* py_U1 = PyList_New(U1.size());
    PyObject* py_U2 = PyList_New(U2.size());

    for (size_t i = 0; i < Z1.size(); ++i) {
        PyList_SetItem(py_Z1, i, PyFloat_FromDouble(Z1[i]));
        PyList_SetItem(py_Z2, i, PyFloat_FromDouble(Z2[i]));
        PyList_SetItem(py_U1, i, PyFloat_FromDouble(U1[i]));
        PyList_SetItem(py_U2, i, PyFloat_FromDouble(U2[i]));
	}
	
	//Import command from Python and run
	const char* command = R"(
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

def plot_data(Z1, Z2, U1, U2):
    # Scatter plot for Standard Normal Distribution
    plt.figure(1)
    plt.scatter(Z1, Z2, color='blue', marker='o', s=1)
    plt.xlabel('Z1')
    plt.ylabel('Z2')
    plt.title('Scatter plot for Standard Normal Distribution')
    plt.show(block=False)  
    
    # Scatter plot for Standard Uniform Distribution
    plt.figure(2)
    plt.scatter(U1, U2, color='red', marker='o', s=1)
    plt.xlabel('U1')
    plt.ylabel('U2')
    plt.title('Scatter plot for Standard Uniform Distribution')
    plt.show(block=False)  
    
    # Histogram plot for Z1
    plt.figure(3)
    plt.hist(Z1, density=True, bins=100, color='skyblue', edgecolor='black')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    plt.plot(x, norm.pdf(x, 0, 1), 'k', linewidth=2)
    plt.xlabel('Z1')
    plt.ylabel('Numbers')
    plt.title('Histogram plot for Z1')
    plt.show(block=False)  
    
    # Histogram plot for Z2
    plt.figure(4)
    plt.hist(Z2, density=True, bins=100, color='skyblue', edgecolor='black')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    plt.plot(x, norm.pdf(x, 0, 1), 'k', linewidth=2)
    plt.xlabel('Z2')
    plt.ylabel('Numbers')
    plt.title('Histogram plot for Z2')
    plt.show(block=False) 

    input("Press Enter to close all plots...")
)";

PyRun_SimpleString(command);


    
    
    
    PyObject* plot_data_func = PyObject_GetAttrString(PyImport_AddModule("__main__"), "plot_data");

    if (plot_data_func && PyCallable_Check(plot_data_func)) {
        // Call the Python function with the arguments
        PyObject* pArgs = PyTuple_Pack(4, py_Z1, py_Z2, py_U1, py_U2);
        PyObject_CallObject(plot_data_func, pArgs);
        Py_DECREF(pArgs);
    } else {
        std::cerr << "Cannot find function 'plot_data'" << std::endl;
    }

    Py_XDECREF(plot_data_func);
    Py_DECREF(py_Z1);
    Py_DECREF(py_Z2);
}


int main() {
	//Importing data from CSV
	string filename_StdNormal = "Data from standard normal distribution.csv";
	string filename_StdUniform = "Data from standard uniform distribution.csv";
	
	vector<double> Z1, Z2; //Standard normal distribution
	vector<double> U1, U2; //standard uniform distribution
	ReadCSV(filename_StdNormal, Z1, Z2);
	ReadCSV(filename_StdUniform, U1, U2);
	
	
	//Python code for plotting
    Py_Initialize();
    call_Python(Z1, Z2, U1, U2);
    Py_Finalize();

    return 0;
    
    
}
