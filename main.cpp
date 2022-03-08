#include <iostream>
#include "dataAcquisition.h"
#include "RTree.h"

int main() {
    RTree rTree;

    rTree.insert(readFile("./../dataTrackMatch.csv"));

    return 0;
}
