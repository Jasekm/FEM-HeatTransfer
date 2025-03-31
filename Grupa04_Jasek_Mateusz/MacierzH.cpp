#include "MacierzH.h"
MacierzH::MacierzH() {};
MacierzH::MacierzH(const ElemUniv& elem, const Jakobian& jacobian, int npc, const wezly& gauss, GlobalData& data,const grid& grid,int elemNumber) {
    macierz.resize(4, vector<double>(4));
    macierzC.resize(4, vector<double>(4));
    macierzHbc.resize(4, vector<double>(4));
    dN_dx.resize(npc, vector<double>(4));
    dN_dy.resize(npc, vector<double>(4));
    P.resize(4);
    for (int i = 0; i < npc; i++) {
        vector<double> temp_dN_dx(4);
        vector<double> temp_dN_dy(4);

        dN_dx[i][0] = jacobian.J1[i][0] * elem.dN_dKsi[i][0] + jacobian.J1[i][1] * elem.dN_dEta[i][0];
        dN_dx[i][1] = jacobian.J1[i][0] * elem.dN_dKsi[i][1] + jacobian.J1[i][1] * elem.dN_dEta[i][1];
        dN_dx[i][2] = jacobian.J1[i][0] * elem.dN_dKsi[i][2] + jacobian.J1[i][1] * elem.dN_dEta[i][2];
        dN_dx[i][3] = jacobian.J1[i][0] * elem.dN_dKsi[i][3] + jacobian.J1[i][1] * elem.dN_dEta[i][3];

        temp_dN_dx[0] = dN_dx[i][0];
        temp_dN_dx[1] = dN_dx[i][1];
        temp_dN_dx[2] = dN_dx[i][2];
        temp_dN_dx[3] = dN_dx[i][3];

        dN_dy[i][0] = jacobian.J1[i][2] * elem.dN_dKsi[i][0] + jacobian.J1[i][3] * elem.dN_dEta[i][0];
        dN_dy[i][1] = jacobian.J1[i][2] * elem.dN_dKsi[i][1] + jacobian.J1[i][3] * elem.dN_dEta[i][1];
        dN_dy[i][2] = jacobian.J1[i][2] * elem.dN_dKsi[i][2] + jacobian.J1[i][3] * elem.dN_dEta[i][2];
        dN_dy[i][3] = jacobian.J1[i][2] * elem.dN_dKsi[i][3] + jacobian.J1[i][3] * elem.dN_dEta[i][3];

        temp_dN_dy[0] = dN_dy[i][0];
        temp_dN_dy[1] = dN_dy[i][1];
        temp_dN_dy[2] = dN_dy[i][2];
        temp_dN_dy[3] = dN_dy[i][3];

        vector<vector<double>> transposed_dN_dx = transpose(temp_dN_dx);
        vector<vector<double>> result_x = multiply(temp_dN_dx, transposed_dN_dx);
        vector<vector<double>> transposed_dN_dy = transpose(temp_dN_dy);
        vector<vector<double>> result_y = multiply(temp_dN_dy, transposed_dN_dy);

        vector<vector<double>> result_xy = add(result_x, result_y);
        result_xy = multiplyByScalar(result_xy, data.Conductivity);
        result_xy = multiplyByScalar(result_xy, jacobian.detJ[i]);
       
        int indexOfWeight1 = i % (npc / static_cast<int>(sqrt(npc)));
        int indexOfWeight2 = i / (npc / sqrt(npc));

        double weight1 = gauss.weights[indexOfWeight1];
        double weight2 = gauss.weights[indexOfWeight2];
        double weight = weight1 * weight2;
        
        result_xy = multiplyByScalar(result_xy, weight);
      //  cout << i + 1 << " pc" << endl;
       // print(result_xy);
       
        
        macierz = add(macierz, result_xy);

        
        vector<vector<double>> transposed_N = transpose(elem.N[i]);
        
        vector<vector<double>> result_N = multiply(elem.N[i], transposed_N);
        result_N = multiplyByScalar(result_N, data.SpecificHeat);
        result_N = multiplyByScalar(result_N, data.Density);
        result_N = multiplyByScalar(result_N, jacobian.detJ[i]);
        result_N = multiplyByScalar(result_N, weight);
        
        macierzC = add(macierzC, result_N);
        
        
        
    }
    for (int k = 0; k < 4; k++) {
        int nextNode = (k + 1) % 4;
       // cout << grid.elements[elemNumber].ID[k] << endl;
        if (grid.nodes[grid.elements[elemNumber].ID[k]-1].BC == 1 && grid.nodes[grid.elements[elemNumber].ID[nextNode]-1].BC == 1) {
            //cout << "Id:" << grid.elements[elemNumber].ID[k]  << endl;
           // cout << "Id nastepne:" << grid.elements[elemNumber].ID[nextNode]  << endl;
           
            vector<vector<double>>  result;
            result.resize(4, vector<double>(4));
            vector<vector<double>> result_N;
            result_N.resize(4, vector<double>(4));

            vector<double> resultP;
            resultP.resize(4);
           vector<double> result_fP;
            result_fP.resize(4);

          for (int j = 0; j < sqrt(npc); j++) {

                    vector<vector<double>> transposed_N = transpose(elem.surface[k].N[j]);
                     result_N = multiply(elem.surface[k].N[j], transposed_N);
                    result_N = multiplyByScalar(result_N, data.Alfa);
                    double weight = gauss.weights[j];
                    result_N = multiplyByScalar(result_N, weight);
                    //print(result_N);
                    result = add(result, result_N);

                    for (int i = 0; i < 4; i++) {
                        resultP[i] = elem.surface[k].N[j][i];
                    }
                    
                    resultP = multiplyVectorByScalar(resultP, data.Tot);
                    resultP = multiplyVectorByScalar(resultP, weight);
                    result_fP = addVectors(result_fP,resultP);

          }
          double detJ = (sqrt(pow(grid.nodes[grid.elements[elemNumber].ID[k] - 1].x - grid.nodes[grid.elements[elemNumber].ID[nextNode] - 1].x, 2) + pow(grid.nodes[grid.elements[elemNumber].ID[k] - 1].y - grid.nodes[grid.elements[elemNumber].ID[nextNode] - 1].y, 2))) / 2;
          result = multiplyByScalar(result, detJ); 
         // print(result);
          result_fP = multiplyVectorByScalar(result_fP, detJ);
          result_fP = multiplyVectorByScalar(result_fP, data.Alfa);
          macierzHbc = add(macierzHbc, result);
          P = addVectors(P, result_fP);
        }

    }
    //cout << "Macierz H: " << endl;
    //print(macierz);
 
    macierz = add(macierz, macierzHbc);
    
   /* cout << "Wektor P: " << endl;
    printVector(P);*/
   
}

vector<vector<double>> MacierzH::transpose(const vector<double>& vec) {
    int size = vec.size();
    vector<vector<double>> transposed(size, vector<double>(1));

    for (int i = 0; i < size; i++) {
        transposed[i][0] = vec[i];
    }

    return transposed;
}

vector<vector<double>> MacierzH::multiply(const vector<double>& vec, const vector<vector<double>>& transposed_vec) {
    int size = vec.size();
    vector<vector<double>> result(size, vector<double>(size, 0));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            result[i][j] = vec[i] * transposed_vec[j][0];
        }
    }

    return result;
}

vector<vector<double>> MacierzH::add(const vector<vector<double>>& A, const vector<vector<double>>& B) {
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

vector<vector<double>> MacierzH::multiplyByScalar(const vector<vector<double>>& matrix, double scalar) {
    vector<vector<double>> result(matrix.size(), vector<double>(matrix[0].size(), 0));

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            result[i][j] = matrix[i][j] * scalar;
        }
    }

    return result;
}
vector<double> MacierzH::multiplyVectorByScalar(const std::vector<double>& vec, double number) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * number;
    }
    return result;
}
vector<double> MacierzH::addVectors(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    // Sprawdzenie, czy wektory maj¹ taki sam rozmiar
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Wektory musz¹ mieæ ten sam rozmiar, aby je dodaæ.");
    }

    std::vector<double> result(vec1.size(), 0);

    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }

    return result;
}

void MacierzH::print_dN_dx() const {
    cout << "dN/dx:" << endl;
    for (const auto& row : dN_dx) {
        cout << "[" << row[0] << ", " << row[1] << ", " << row[2] << ", " << row[3] << "]" << endl;
    }
}

void MacierzH::print_dN_dy() const {
    cout << "dN/dy:" << endl;
    for (const auto& row : dN_dy) {
        cout << "[" << row[0] << ", " << row[1] << ", " << row[2] << ", " << row[3] << "]" << endl;
    }
}
void MacierzH::print_macierzH() const {
    cout << "MacierzH" << endl;
    for (const auto& row : macierz) {
        cout << "[" << row[0] << ", " << row[1] << ", " << row[2] << ", " << row[3] << "]" << endl;
    }
}
void MacierzH::print_macierzC() const {
    cout << "MacierzC" << endl;
    for (const auto& row : macierzC) {
        cout << "[" << row[0] << ", " << row[1] << ", " << row[2] << ", " << row[3] << "]" << endl;
    }
}

void MacierzH::print(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void MacierzH::printVector(const std::vector<double>& vec) {
    for (double val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}