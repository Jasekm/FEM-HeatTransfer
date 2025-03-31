#include "GlobalData.h"
#include "shared_data.h"

GlobalData GloadData() {
    ifstream file(FILENAME);

    if (!file.is_open()) {
        cout << "Unable to open file" << endl;
        exit(1);
    }

    string line;
    int SimulationTime = 0, SimulationStepTime = 0, Conductivity = 0, Alfa = 0, Tot = 0,
        InitialTemp = 0, Density = 0, SpecificHeat = 0, nN = 0, nE = 0, Npc = 0;

    while (getline(file, line)) {
        if (line.find("SimulationTime") != string::npos) {
            SimulationTime = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("SimulationStepTime") != string::npos) {
            SimulationStepTime = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("Conductivity") != string::npos) {
            Conductivity = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("Alfa") != string::npos) {
            Alfa = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("Tot") != string::npos) {
            Tot = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("InitialTemp") != string::npos) {
            InitialTemp = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("Density") != string::npos) {
            Density = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("SpecificHeat") != string::npos) {
            SpecificHeat = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("Nodes number") != string::npos) {
            nN = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("Elements number") != string::npos) {
            nE = stoi(line.substr(line.find_last_of(" ") + 1));
        }
        else if (line.find("Npc") != string::npos) {
            Npc = stoi(line.substr(line.find_last_of(" ") + 1));
        }
    }

    file.close();

    return GlobalData(SimulationTime, SimulationStepTime, Conductivity, Alfa, Tot, InitialTemp, Density, SpecificHeat, nN, nE,Npc);
}
void displayData(GlobalData data) {
    cout << "Simulation Time: " << data.SimulationTime << endl;
    cout << "Simulation Step Time: " << data.SimulationStepTime << endl;
    cout << "Conductivity: " << data.Conductivity << endl;
    cout << "Alfa: " << data.Alfa << endl;
    cout << "Tot: " << data.Tot << endl;
    cout << "Initial Temp: " << data.InitialTemp << endl;
    cout << "Density: " << data.Density << endl;
    cout << "Specific Heat: " << data.SpecificHeat << endl;
    cout << "nN (Nodes number): " << data.nN << endl;
    cout << "nE (Elements number): " << data.nE << endl;
    cout << "Npc: " << data.Npc << endl;


}