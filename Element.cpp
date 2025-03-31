#include "Element.h"
#include "shared_data.h"

vector<element> EloadData(int nE) {
    ifstream file(FILENAME);

    if (!file.is_open()) {
        cout << "Unable to open file" << endl;
        exit(1);
    }

    string line;
    vector<element> elements;
    bool foundElementSection = false;
    int elementsCounter = 0;
    while (getline(file, line) && elementsCounter < nE) {
        if (line.find("*Element") != string::npos) {
            foundElementSection = true;
            continue;
        }

        if (foundElementSection) {
            stringstream ss(line);
            char comma;
            int id = 0;
            int a = 0;
            int b = 0;
            int c = 0;
            int d = 0;
            ss >> id >> comma >> a >> comma >> b >> comma >> c >> comma >> d;
            elements.emplace_back(a, b, c, d);
            elementsCounter++;
        }
    }

    file.close();
    return elements;
}
