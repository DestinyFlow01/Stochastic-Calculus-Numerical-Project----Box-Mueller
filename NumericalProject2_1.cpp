#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<vector>
#include<algorithm>
#include<fstream>
#include<cmath>
#include <Python.h>
#include <sstream>
#include <string>
using namespace std;

//Details in generating uniform distribution
int low = 0, high = 50000000, N = 10000;
double range_low = 1e-8, range_high = 2;
	
	
	
	
const double PI = 3.14159265;	

//For generate random number with uniform distribution (only for generating)
class RandomNumber {
	private : 
		vector<double> Random_Uniform;
		vector<double> Random_StdUniform;
		int low, high, N;	
		double range_low, range_high;	
	public :
		RandomNumber(int low, int high, int N, double range_low, double range_high); //Constructor (Generating N random numbers with uniform distribution, generated from integer)
		void Print_Uniform();
		void Print_StdUniform();
		vector<double> Access_RandomUniform();
		vector<double> Access_RandomStdUniform();
};

class RandomVector {
	private : 
		int dim;
		vector<RandomNumber> vec_Random_Uniform;
		
		vector<vector<double>> vec_StdNormal;
	public : 
		//Generating and Printing Random Vector with uniform distribution
		RandomVector(int dim);
		void Print_vecUniform();
		void Print_vecStdUniform();
		
		//Converting from standard normal to standard Gaussian using Box-Muller method
		void StdUniform2StdNormal();
		void Print_vecStdNormal();
		
		//Output in CSV
		void OutputCSV_Uniform(string name);
		void OutputCSV_StdUniform(string name);
		void OutputCSV_StdNormal(string name);
		
		//Finding covariance
		double Covariance_StdNormal();
};


//Defining member function for RandomNumber
RandomNumber::RandomNumber(int low, int high, int N, double r_low, double r_high) : low(low), high(high), N(N), range_low(r_low), range_high(r_high) {
	
	
	int temp, value;
	vector<double> Random_Number;
	//Generating random number from integer
	for(int i = 0; i<N; i++) {
		temp = rand();
		
		//Converting to the wanted range
		value = temp%(high-low+1) + low;
		
		//Inserting it to Random_Uniform first
		Random_Number.push_back(double(value));
	}
	
	
	//Finding maximum value of Random_Uniform in the range
	auto max_Element = max_element(Random_Number.begin(), Random_Number.end());
	
	Random_Uniform.resize(Random_Number.size());
	Random_StdUniform.resize(Random_Number.size());
	
	//Dividing them into 0 and 1, and change the range from A to B
	transform(Random_Number.begin(), Random_Number.end(), Random_Uniform.begin(), [max_Element, r_low, r_high](double x) {return (x/ *max_Element)*(r_high - r_low) + r_low;});
	transform(Random_Uniform.begin(), Random_Uniform.end(), Random_StdUniform.begin(), [r_low, r_high](double x) {return x/(r_high - r_low); });
	
	max_Element = max_element(Random_StdUniform.begin(), Random_StdUniform.end());
	auto min_Element = min_element(Random_StdUniform.begin(), Random_StdUniform.end());
	cout<<"Random number from uniform distribution has been generated with min element = "<<*min_Element<<" and max element = "<<*max_Element<<endl;
}

void RandomNumber::Print_Uniform() {
	//Outputting first 10 and last 10 if N>20. if N<20, output all
	cout<<"Printing samples from uniform distribution \n";
	if(Random_Uniform.size()>20) {
		cout<<"First 10 : "<<setprecision(6);
		for(int i = 0; i<10; i++) {
			cout<<Random_Uniform[i]<<", ";
		}
		
		cout<<"\nLast 10 : ";
		for(int i = Random_Uniform.size()-10; i<Random_Uniform.size(); i++) {
			cout<<Random_Uniform[i]<<", ";
		}
		cout<<endl;
	}
	else {
		cout<<"Uniform distribution sample : "<<setprecision(6);
		for(int i = 0; i<N; i++) {
			cout<<Random_Uniform[i]<<", ";
		}
	}
	cout<<endl;
}

void RandomNumber::Print_StdUniform() {
	//Outputting first 10 and last 10 if N>20. if N<20, output all
	cout<<"Printing samples from standard uniform distribution \n";
	if(Random_StdUniform.size()>20) {
		cout<<"First 10 : "<<setprecision(6);
		for(int i = 0; i<10; i++) {
			cout<<Random_StdUniform[i]<<", ";
		}
		
		cout<<"\nLast 10 : ";
		for(int i = Random_StdUniform.size()-10; i<Random_StdUniform.size(); i++) {
			cout<<Random_StdUniform[i]<<", ";
		}
		cout<<endl;
	}
	else {
		cout<<"Uniform distribution sample : "<<setprecision(6);
		for(int i = 0; i<N; i++) {
			cout<<Random_StdUniform[i]<<", ";
		}
	}
	cout<<endl;
}

vector<double> RandomNumber::Access_RandomStdUniform() { return Random_StdUniform; }
vector<double> RandomNumber::Access_RandomUniform() { return Random_Uniform; }

//Defining member function for RandomVector
RandomVector::RandomVector(int dim) : dim(dim) {
	cout<<"Generating Random Vector from uniform distribution : \n";
	//Generating Random Number and combine into vector
	for(int i = 0; i<dim; i++) {
		RandomNumber U(low, high, N, range_low, range_high);
		vec_Random_Uniform.push_back(U);
	}
	
	
}

void RandomVector::Print_vecUniform() {
	cout<<"\nUniform distribution : \n";
	for(int i = 0; i<dim; i++) {
		cout<<"U"<<i+1<<" : \n";
		vec_Random_Uniform[i].Print_Uniform();
		cout<<endl;
	}
}

void RandomVector::Print_vecStdUniform() {
	cout<<"\nStandard Uniform distribution : \n";
	for(int i = 0; i<dim; i++) {
		cout<<"U"<<i+1<<" : \n";
		vec_Random_Uniform[i].Print_StdUniform();
		cout<<endl;
	}

}

void RandomVector::StdUniform2StdNormal() {
	vector<double> Z1, Z2;
	vector<double> U1 = vec_Random_Uniform[0].Access_RandomStdUniform();
	vector<double> U2 = vec_Random_Uniform[1].Access_RandomStdUniform();
	
	for(int i = 0; i<N; i++) {	
		double z1, z2;
		z1 = sqrt(-2*log(U1[i]))*cos(2*PI*U2[i]);
		z2 = sqrt(-2*log(U1[i]))*sin(2*PI*U2[i]);
		
		if(isnan(z1) or isnan(z2)) {
			Z1.push_back(0); Z2.push_back(0);
		}
		else {
			Z1.push_back(z1); Z2.push_back(z2);
		}
		
	}
	
	vec_StdNormal.push_back(Z1); vec_StdNormal.push_back(Z2); 
}

void RandomVector::Print_vecStdNormal() {
	cout<<"\nStandard Normal distribution : \n";
	
	//Outputting first 10 and last 10 if N>20. if N<20, output all
	cout<<"Printing Converted from standard normal distribution \n";
	for(int j = 0; j<dim; j++) {
		if(vec_StdNormal[j].size()>20) {
			cout<<"First 10 : "<<setprecision(6);
			for(int i = 0; i<10; i++) {
				cout<<vec_StdNormal[j][i]<<", ";
			}
			
			cout<<"\nLast 10 : ";
			for(int i = vec_StdNormal[j].size()-10; i<vec_StdNormal[j].size(); i++) {
			cout<<vec_StdNormal[j][i]<<", ";
			}
			cout<<endl;
		}
		else {
			cout<<"Standard Normal distribution sample : "<<setprecision(6);
			for(int i = 0; i<vec_StdNormal[j].size(); i++) {
				cout<<vec_StdNormal[j][i]<<", ";
			}
		}
	}
	
	cout<<endl;
}

void RandomVector::OutputCSV_Uniform(string name) {
	vector<double> U1 = vec_Random_Uniform[0].Access_RandomUniform(), U2 = vec_Random_Uniform[1].Access_RandomUniform();
	
	ofstream ofs;
	ofs.open(name);
	ofs << "Data,U1,U2\n";
	for (int i = 0; i<U1.size(); i++) {
		ofs << i << "," << U1[i]<< "," << U2[i] << "\n";
	}
	ofs.close();
} 

void RandomVector::OutputCSV_StdUniform(string name) {
	vector<double> U1 = vec_Random_Uniform[0].Access_RandomStdUniform(), U2 = vec_Random_Uniform[1].Access_RandomStdUniform();
	
	ofstream ofs;
	ofs.open(name);
	ofs << "Data,U1,U2\n";
	for (int i = 0; i<U1.size(); i++) {
		ofs << i << "," << U1[i]<< "," << U2[i] << "\n";
	}
	ofs.close();
} 

void RandomVector::OutputCSV_StdNormal(string name) {
	ofstream ofs;
	ofs.open(name);
	ofs << "Data,Z1,Z2\n";
	for (int i = 0; i<N; i++) {
		ofs << i << "," << vec_StdNormal[0][i]<< "," << vec_StdNormal[1][i] << "\n";
	}
	ofs.close();
} 

//Finding covariance
double RandomVector::Covariance_StdNormal() {
	double sum = 0;
	
	for(int i = 0; i<N; i++) {
		sum += vec_StdNormal[0][i]*vec_StdNormal[1][i];
	}
	
	double covariance = sum/N;
	return covariance;
}


//Plotting 
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
    plt.hist(Z1, density=True, bins=100, color='skyblue', edgecolor='black', label='Z1')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    plt.plot(x, norm.pdf(x, 0, 1), 'k', linewidth=2, label = 'Normal N(0,1)')
    plt.xlabel('Z1')
    plt.ylabel('Numbers')
    plt.title('Histogram plot for Z1')
    plt.legend()
    plt.show(block=False)  
    
    # Histogram plot for Z2
    plt.figure(4)
    plt.hist(Z2, density=True, bins=100, color='skyblue', edgecolor='black', label='Z2')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    plt.plot(x, norm.pdf(x, 0, 1), 'k', linewidth=2, label = 'Normal N(0,1)')
    plt.xlabel('Z2')
    plt.ylabel('Numbers')
    plt.legend()
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




//Main function
int main() {
	srand(time(0)); //initializing random seed
	cout<<"Random Vector U = (U1, U2) : \n";
	RandomVector RV1(2);
//	RV1.Print_vecUniform();
//	RV1.Print_vecStdUniform();
	
	//Convert to Normal distribution using Box-Mueller method
	cout<<"Converting from Standard Uniform to Standard Normal distribution : \n";
	RV1.StdUniform2StdNormal();
	cout<<"Converting done\n\n";
	
	//Printing result : 
//	RV1.Print_vecStdNormal();
	
	//Outputting result
	RV1.OutputCSV_Uniform("Data from uniform distribution.csv");
	cout<<"Data from uniform distribution have been outputted\n";
	
	RV1.OutputCSV_StdUniform("Data from standard uniform distribution.csv");
	cout<<"Data from standard uniform distribution have been outputted\n";
	
	RV1.OutputCSV_StdNormal("Data from standard normal distribution.csv");
	cout<<"Data from standard normal distribution have been outputted\n";
	
	//Finding covariance
	double covariance = RV1.Covariance_StdNormal();
	cout<<"\nCovariance of Z1 and Z2 = "<<covariance<<endl;
	
	
	//Plotting process
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
    
    /*
    Type this in CMD
    g++ NumericalProject2_1.cpp -IC:\Python313\include -LC:\Python313\libs -lpython313 -o NumericalProject2_1 -Wl,--enable-auto-import
	NumericalProject2_1.exe
    */
}