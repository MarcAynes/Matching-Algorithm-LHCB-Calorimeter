#ifndef MATCHINGALGORITHM_RTREE_H
#define MATCHINGALGORITHM_RTREE_H

#include "hit.h"
#include <vector>

using namespace std;

/*
 * R-Tree Data structure
 * spatial indexing of data
 */

typedef struct Node{
    void* node;
    double xMax;
    double yMax;
    double xMin;
    double yMin;
    bool leaf;              //bool leaf it's used to know if the next node (going down) it's a leaf node or a regular node
}Node;

typedef struct leafNode{
    hit Hit;
    double X;
    double Y;
}leafNode;

class RTree {
private:


public:
    /*
     * Insert Function
     * add all data to the data structure,
     * Using recursive function be aware of the stack used and the stack available
     */

    void insert(vector<vector<hit>> data){

    }
};


#endif //MATCHINGALGORITHM_RTREE_H
