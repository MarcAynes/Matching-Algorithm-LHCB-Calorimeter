#include <iostream>
#include "dataAcquisition.h"
#include "RTree.h"
#include "trace.h"

int main() {
    RTree rTree;
    vector<vector<trace>> traces;

    rTree.insert(readFile("./../dataTrackMatch.csv", &traces));

    vector<vector<vector<hit>>> possibleHits = rTree.searchWithTrace(traces);

    return 0;
}
