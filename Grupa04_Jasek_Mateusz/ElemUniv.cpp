#include "ElemUniv.h"
#include <iostream>

using namespace std;  

ElemUniv::ElemUniv(int npc, const wezly& gauss) {
    dN_dKsi.resize(npc, vector<double>(4));
    dN_dEta.resize(npc, vector<double>(4));
    N.resize(npc, vector<double>(4));
    surface.resize(4);
    for (int i = 0; i < 4; i++) {
        surface[i].N.resize(sqrt(npc), vector<double>(4));
    }
    for (int i = 0; i < npc; i++) {
        int nodeValueEta = i / (npc / sqrt(npc));
        int nodeValueKsi = i % (npc / static_cast<int>(sqrt(npc)));

        for (int j = 0; j < 4; j++) {
            if (j == 0) {
                dN_dKsi[i][j] = -0.25 * (1 - gauss.nodes[nodeValueEta]);
                dN_dEta[i][j] = -0.25 * (1 - gauss.nodes[nodeValueKsi]);
                N[i][j] = 0.25 * (1 - gauss.nodes[nodeValueKsi]) * (1 - gauss.nodes[nodeValueEta]);
            }
            else if (j == 1) {
                dN_dKsi[i][j] = 0.25 * (1 - gauss.nodes[nodeValueEta]);
                dN_dEta[i][j] = -0.25 * (1 + gauss.nodes[nodeValueKsi]);
                N[i][j] = 0.25 * (1 + gauss.nodes[nodeValueKsi]) * (1 - gauss.nodes[nodeValueEta]);

            }
            else if (j == 2) {
                dN_dKsi[i][j] = 0.25 * (1 + gauss.nodes[nodeValueEta]);
                dN_dEta[i][j] = 0.25 * (1 + gauss.nodes[nodeValueKsi]);
                N[i][j] = 0.25 * (1 + gauss.nodes[nodeValueKsi]) * (1 + gauss.nodes[nodeValueEta]);

            }
            else if (j == 3) {
                dN_dKsi[i][j] = -0.25 * (1 + gauss.nodes[nodeValueEta]);
                dN_dEta[i][j] = 0.25 * (1 - gauss.nodes[nodeValueKsi]);
                N[i][j] = 0.25 * (1 - gauss.nodes[nodeValueKsi]) * (1 + gauss.nodes[nodeValueEta]);

            }
        }
    }
    for (int k = 0; k < 4; k++) {

        for (int i = 0; i < sqrt(npc); i++) {
            int nodeValueEta = i;
            int nodeValueKsi = i;
            double ksi = 0;
            double eta = 0;
            if (k == 0) {
                 ksi = gauss.nodes[nodeValueKsi];
                 eta = -1;
            }
            else if (k == 1) {
                 ksi = 1;
                 eta = gauss.nodes[nodeValueEta];
            }
            else if (k == 2) {
                 ksi = gauss.nodes[nodeValueKsi];
                 eta = 1;
            }
            else if (k == 3) {
                 ksi = -1;
                 eta = gauss.nodes[nodeValueEta];
            }
            
            
            for (int j = 0; j < 4; j++) {
                if (j == 0) {
                    surface[k].N[i][j] = 0.25 * (1 - ksi) * (1 - eta);

                }
                else if (j == 1) {
                    surface[k].N[i][j] = 0.25 * (1 + ksi) * (1 - eta);

                }
                else if (j == 2) {
                    surface[k].N[i][j] = 0.25 * (1 + ksi) * (1 + eta);

                }
                else if (j == 3) {
                    surface[k].N[i][j] = 0.25 * (1 - ksi) * (1+ eta);

                }
            }
        }
    }
}

void ElemUniv::print() {
    cout << "dN_dKsi:" << endl;
    for (const auto& row : dN_dKsi) {
        for (const auto& val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    cout << "dN_dEta:" << endl;
    for (const auto& row : dN_dEta) {
        for (const auto& val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    for (int k = 0; k < 4; ++k) {
        cout << "Surface " << k + 1 << ":" << endl;
        for (const auto& row : surface[k].N) {
            for (const auto& val : row) {
                cout << val << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << "N:" << endl;
    for (const auto& row : N) {
        for (const auto& val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}
