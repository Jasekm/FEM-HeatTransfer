#include "Node.h"
#include "shared_data.h"

vector<int> LoadBC() {
    ifstream file(FILENAME);

    if (!file.is_open()) {
        cout << "Unable to open file" << endl;
        exit(1);
    }

    string line;
    vector<int> BC;
    bool foundBCSection = false;

    while (getline(file, line)) {
        if (line.find("*BC") != string::npos) {
            foundBCSection = true;
            continue;
        }

        if (foundBCSection) {
            stringstream ss(line);
            int node;
            while (ss >> node) {
                BC.push_back(node);
                if (ss.peek() == ',') ss.ignore();  
            }
            break;  
        }
    }

    file.close();
    return BC;
}

vector<node> NloadData(int nN) {
    ifstream file(FILENAME);

    if (!file.is_open()) {
        cout << "Unable to open file" << endl;
        exit(1);
    }

    string line;
    vector<int> BC = LoadBC();  
    vector<node> nodes;
    bool foundNodeSection = false;
    int nodesCounter = 0;

    while (getline(file, line) && nodesCounter < nN) {
        if (line.find("*Node") != string::npos) {
            foundNodeSection = true;
            continue;
        }

        if (foundNodeSection) {
            stringstream ss(line);
            char comma;
            double x = 0;
            double y = 0;
            int id = 0;
            ss >> id >> comma >> x >> comma >> y;

            int bc = (find(BC.begin(), BC.end(), id) != BC.end()) ? 1 : 0;

            nodes.emplace_back(x, y, bc);
            nodesCounter++;
        }
    }

    file.close();
    return nodes;
}
