#include "Grid.h"
#include <iomanip>

void displayGrid(grid newGrid) {
    int counterN = 1;
    cout << "Nodes: " << endl;
    for (const node& node : newGrid.nodes) {
       // cout << setprecision(11) << counterN << " x = " << node.x << ", y = " << node.y << endl;
        cout << setprecision(11)  << " x = " << node.x << ", y = " << node.y << endl;

        counterN++;
    }
    int counterE = 1;
    cout << "Elements: " << endl;
    for (const element& element : newGrid.elements) {
      //  cout << counterE << ", " << element.ID[0] << ", " << element.ID[1] << ", " << element.ID[2] << ", " << element.ID[3] << endl;
        cout   << element.ID[0] << ", " << element.ID[1] << ", " << element.ID[2] << ", " << element.ID[3] << endl;

        counterE++;
    }
}
