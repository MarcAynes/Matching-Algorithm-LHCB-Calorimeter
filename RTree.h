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
    short numNode;          //it could be a char for more memory optimization but the way C structures work it won't change the memory allocated
    bool leaf;              //bool leaf it's used to know if the next node (going down) it's a leaf node or a regular node
}Node;

typedef struct leafNode{
    hit Hit;
    double X;
    double Y;
}leafNode;

class RTree {
private:
    Node root{};

    /*
     * Recursive function where we will insert the leaf Node
     */

    void insert(leafNode leaf, Node* current){
        if (current->leaf){                 //root is leaf
            if (current->numNode < 5){      //if there is space
                if (current->numNode == 0){ //it does not have info insered
                    current->node = (leafNode*) malloc(sizeof (leafNode));
                    current->numNode++;
                    *((leafNode*) (current->node)) = leaf;
                }
            }
        }
    }

public:
    /*
     * Insert Function
     * add all data to the data structure,
     * Using recursive function be aware of the stack used and the stack available
     */

    void insert(vector<vector<hit>> data){
        for (int i = 0; i < data.size(); i++){
            for (int j = 0; j < data[i].size(); j++) {
                leafNode newLeaf;
                newLeaf.X = data[i][j].X;
                newLeaf.Y = data[i][j].Y;
                newLeaf.Hit = data[i][j];
                insert(newLeaf, &root);
            }
        }
    }

    RTree(){
        root.leaf = true;
        root.node = nullptr;
        root.numNode = 0;
        root.xMax = -1;
        root.xMin = -1;
        root.yMax = -1;
        root.yMin = -1;
    }
};


#endif //MATCHINGALGORITHM_RTREE_H
