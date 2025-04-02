#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include "GlobalData.h"
#include "Node.h"
#include "Element.h"
#include "Grid.h"
#include "Wezly.h"
#include "ElemUniv.h"
#include "Jakobian.h"
#include "MacierzH.h"



using namespace std;

double funkcja_kwadratowa(double x) {
    return 5 * x * x + 3 * x + 6;
}
double funkcja_kwadratowa2(double x, double y) {
    return 5 * x * x * y * y + 3 * x * y + 6;
}
void metoda_kwadratur1d(double (*funkcja)(double), const wezly& gauss) {
    double suma = 0;
    for (int i = 0; i <= gauss.DEGREE; i++) {
        suma += gauss.weights[i] * funkcja(gauss.nodes[i]);
    }
    cout << "Dla N: " << gauss.DEGREE << " = " << suma << endl;
}
void metoda_kwadratur2d(double (*funkcja)(double, double), const wezly& gauss) {
    double suma = 0;
    for (int i = 0; i <= gauss.DEGREE; i++) {
        for (int j = 0; j <= gauss.DEGREE; j++) {
            suma += gauss.weights[i] * gauss.weights[j] * funkcja(gauss.nodes[i], gauss.nodes[j]);
        }
    }
    cout << "Dla N: " << gauss.DEGREE << " = " << suma << endl;
}
vector<vector<double>> add(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        throw invalid_argument("Macierze musz¹ mieæ te same wymiary.");
    }

    vector<vector<double>> result(A.size(), vector<double>(A[0].size(), 0));

    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}
vector<vector<double>> multiplyByScalar(const vector<vector<double>>& matrix, double scalar) {
    vector<vector<double>> result(matrix.size(), vector<double>(matrix[0].size(), 0));

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            result[i][j] = matrix[i][j] * scalar;
        }
    }

    return result;
}
vector<vector<double>> divideByScalar(const vector<vector<double>>& matrix, double scalar) {
    if (scalar == 0) {
        throw invalid_argument("Nie mo¿na dzieliæ przez zero.");
    }

    vector<vector<double>> result(matrix.size(), vector<double>(matrix[0].size(), 0));

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            result[i][j] = matrix[i][j] / scalar;
        }
    }

    return result;
}
vector<double> multiplyMatrixByVector(const vector<vector<double>>& A, const vector<double>& v) {
    if (A[0].size() != v.size()) {
        throw invalid_argument("Liczba kolumn w macierzy musi byæ równa liczbie elementów w wektorze.");
    }

    size_t m = A.size();  
    size_t n = A[0].size();  

    vector<double> result(m, 0); 

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i] += A[i][j] * v[j];
        }
    }

    return result;
}
vector<vector<double>> addVectorToMatrix(const vector<vector<double>>& A, const vector<double>& v) {
    if (A[0].size() != v.size()) {
        throw invalid_argument("Liczba kolumn w macierzy musi byæ równa liczbie elementów w wektorze.");
    }

    vector<vector<double>> result = A;  

    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            result[i][j] += v[j];  
        }
    }

    return result;
}
vector<double> addMatrixToVector(const vector<vector<double>>& matrix, const vector<double>& vec) {
    if (matrix[0].size() != vec.size()) {
        throw std::invalid_argument("Wektor musi mieæ tyle samo elementów co liczba kolumn macierzy.");
    }

    vector<double> result(vec.size(), 0.0);

    for (size_t j = 0; j < matrix[0].size(); j++) { 
        result[j] = vec[j]; 
        for (size_t i = 0; i < matrix.size(); i++) { 
            result[j] += matrix[i][j];
        }
    }

    return result;
}
vector<double> addVectors(const vector<double>& vec1, const vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Wektory musz¹ mieæ tê sam¹ d³ugoœæ.");
    }

    vector<double> result(vec1.size(), 0.0);
    for (size_t i = 0; i < vec1.size(); i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}


class H_global {
public:
    vector<vector<double>> macierz;
    vector<vector<double>> macierzC;
    
    H_global(const grid& newGrid,int nE, int nN ) {
        macierz.resize(nN, vector<double>(nN));
        macierzC.resize(nN, vector<double>(nN));
        for (int i = 0; i < nE; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    macierz[newGrid.elements[i].ID[j]-1][newGrid.elements[i].ID[k]-1] += newGrid.elements[i].macierzH->macierz[j][k];
                    macierzC[newGrid.elements[i].ID[j]-1][newGrid.elements[i].ID[k]-1] += newGrid.elements[i].macierzH->macierzC[j][k];
                }
            }
        }
        
    }
    void print_MacierzH() const {
        for (const auto& row : macierz) { 
            for (const auto& value : row) { 
                cout << setw(10) << value << " "; 
            }
            cout << endl;
        }
    }
    void print_MacierzC() const {
        for (const auto& row : macierzC) {
            for (const auto& value : row) {
                cout << setw(10) << value << " ";
            }
            cout << endl;
        }
    }
};
class P_global {
public:
    vector<double> vectorP;

    P_global(const grid& newGrid, int nE, int nN) {
        vectorP.resize(nN);
        for (int i = 0; i < nE; i++) {
                for (int k = 0; k < 4; k++) {
                    vectorP[newGrid.elements[i].ID[k] - 1] += newGrid.elements[i].macierzH->P[k];
                }
        }
        
    }
    void print() const {
        std::cout << "vectorP: ";
        for (const double& val : vectorP) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    
};
class global_system_equation {
public:
    vector<vector<double>> macierzH_global;  
    vector<vector<double>> macierzC_global;  
    vector<double> vectorP_global;           
    
    vector<double> vectorResult;

    vector<vector<double>> macierzTemp;
    vector<double> vectorTemp;
    double t = 0;
    vector<double> t0;
    
    global_system_equation(const vector<vector<double>>& macierzH, const vector<vector<double>>& macierzC, const vector<double>& vectorP, GlobalData& data, const grid& grid) {
        macierzH_global = macierzH;
        macierzC_global = macierzC;
        vectorP_global = vectorP;
        vectorResult.resize(vectorP.size());
        macierzTemp = macierzH;
        vectorTemp = vectorP;
        t = data.SimulationStepTime;
        t0.resize(data.nN);
        for (int i = 0; i < data.nN; i++) {
            t0[i] = grid.nodes[i].Temp;
        }
        
    }
    
    void postepowanie_proste() {
        
        for (int k = 0; k < macierzTemp.size() - 1; k++) {
            for (int i = k + 1; i < macierzTemp.size(); i++) {
                float m = macierzTemp[i][k] / macierzTemp[k][k];
                for (int j = 0; j < macierzTemp.size(); j++) {
                    macierzTemp[i][j] -= macierzTemp[k][j] * m;
                }
                vectorTemp[i] -= vectorTemp[k] * m;  
            }
        }

    }

    void postepowanie_odwrotne() {
        
        for (int i = macierzTemp.size() - 1; i >= 0; i--) {
            float x = vectorTemp[i] / macierzTemp[i][i];  
            vectorResult[i] = x;

            for (int j = i - 1; j >= 0; j--) {
                vectorTemp[j] -= macierzTemp[j][i] * x;  
                macierzTemp[j][i] = 0;  
            }
        }
    }
    vector<double> oblicz_uklad() {
        macierzC_global = divideByScalar(macierzC_global, t);
         macierzTemp = add(macierzH_global, macierzC_global);
        vector<double> resultR = multiplyMatrixByVector(macierzC_global, t0);
        vectorTemp = addVectors(resultR, vectorP_global);
        //printVector(vectorTemp);
        postepowanie_proste();
        postepowanie_odwrotne();

        return vectorResult;
    }
   
    void printResult() const {
        std::cout << "Result: ";
        for (const double& val : vectorResult) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    void printMatrix(const vector<vector<double>>& matrix) const {
        std::cout << "Matrix "  << std::endl;
        for (const auto& row : matrix) {
            for (const double& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void printVector(const vector<double>& vec) const {
        std::cout << "Vector " <<endl;
        for (const double& val : vec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
};
int main() {
    GlobalData data = GloadData();

    vector<node> nodes = NloadData(data.nN);
    vector<element> elements = EloadData(data.nE);
    for (int i = 0; i < data.nN; i++) {
        nodes[i].Temp = data.InitialTemp
    }
    grid newGrid(data.nN, data.nE, elements, nodes);
  

    wezly gauss0 = wezly(0);
    wezly gauss1 = wezly(1);
    wezly gauss2 = wezly(2);
    wezly gauss3 = wezly(3);
    wezly gauss4 = wezly(4);

    wezly gauss = wezly(0);

    switch (data.Npc) {
    case 1:
        break;
    case 4:
        gauss = gauss1;
        break;
    case 9:
        gauss = gauss2;
        break;
    case 16:
        gauss = gauss3;
        break;
    default:
        break;
    }


    ElemUniv elem(data.Npc, gauss);
    int iteracje = data.SimulationTime / data.SimulationStepTime;
    for (int j = 0; j < iteracje; j++) {
        for (int i = 0; i < elements.size(); i++) {
            Jakobian jacobian(elem, newGrid, i, data.Npc);
            newGrid.elements[i].jacobian = jacobian;
            MacierzH* macierzH = new MacierzH(elem, jacobian, data.Npc, gauss, data, newGrid, i);
            newGrid.elements[i].macierzH = macierzH;
       
        }
  
        H_global globalH(newGrid, data.nE, data.nN);
        P_global globalP(newGrid, data.nE, data.nN);
       
        global_system_equation system(globalH.macierz, globalH.macierzC, globalP.vectorP, data, newGrid);
        

        vector<double> t1= system.oblicz_uklad();

        for (int i = 0; i < data.nN; i++) {
            newGrid.nodes[i].Temp = t1[i];
            
        }
        double min = t1[0];
        double max = t1[0];
        for (int i = 0; i < data.nN; i++) {
            if (t1[i] < min) {
                min = t1[i];
            }
            else if (t1[i] > max) {
                max = t1[i];
            }
        }
        cout << fixed << setprecision(10) << min << "  " << max << endl;

    }

    return 0;
}
