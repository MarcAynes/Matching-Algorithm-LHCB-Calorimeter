#ifndef MATCHINGALGORITHM_RTREE_H
#define MATCHINGALGORITHM_RTREE_H

#define m 3
#define M 5

#include "hit.h"
#include <vector>
#include <float.h>
#include <bits/stdc++.h>

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

    double euclideanDistance(double x1, double y1, double x2, double y2)
    {
        // Calculating distance
        return sqrt(pow(x2 - x1, 2) +
                    pow(y2 - y1, 2));
    }

    /*
     * Recursive function where we will insert the leaf Node
     */

    void insert(leafNode leaf, Node* current){
        if (current->leaf){                 //root is leaf
            if (current->numNode < M){      //if there is space
                if (current->numNode == 0){ //it does not have info insered
                    current->node = (leafNode*) malloc(sizeof (leafNode));
                    current->numNode++;
                    *((leafNode*) (current->node)) = leaf;
                    current->yMin = leaf.Y; //adding X and Y data
                    current->yMax = leaf.Y;
                    current->xMin = leaf.X;
                    current->xMax = leaf.X;
                }else {                     //does have info inserted
                    current->numNode++;
                    current->node = (leafNode *) realloc(current->node, sizeof(leafNode) * (current->numNode));
                    ((leafNode *) (current->node))[(current->numNode) - 1] = leaf;

                    if (current->xMax < leaf.X) {   //modifying X and Y data if needed
                        current->xMax = leaf.X;
                        return;
                    }

                    if (current->xMin > leaf.X) {
                        current->xMin = leaf.X;
                        return;
                    }

                    if (current->yMax < leaf.Y) {
                        current->yMax = leaf.Y;
                        return;
                    }

                    if (current->yMin > leaf.Y) {
                        current->yMin = leaf.Y;
                        return;
                    }

                    /*                                                  some prints for debug
                    if (current->numNode >= 3){
                        int a = 0;
                        while (a < current->numNode){
                            printf("id - %d\n", ((leafNode*)current->node)[a].Hit.id);
                            a++;
                        }
                    }
                     */
                }
            }else{      // there is not space
                int index1 = 0, index2 = 0;
                double dist = DBL_MAX;
                for (int i = 0; i < current->numNode; i++) {     //O(N^2) cost at finding 2 furthest points, but for 6 points it's not critical
                    for (int j = i + 1; j < current->numNode; j++) {
                        double auxDist = euclideanDistance(((leafNode *) current->node)[i].X,
                                                           ((leafNode *) current->node)[i].Y,
                                                           ((leafNode *) current->node)[j].X,
                                                           ((leafNode *) current->node)[j].Y);
                        if (auxDist < dist) {
                            index1 = i;
                            index2 = j;
                            dist = auxDist;
                        }
                    }
                }

                /*current->numNode++;
                current->node = (leafNode *) realloc(current->node, sizeof(leafNode) * (current->numNode));
                ((leafNode *) (current->node))[(current->numNode) - 1] = leaf;

                if (current->xMax < leaf.X) {   //modifying X and Y data if needed
                    current->xMax = leaf.X;
                    return;
                }

                if (current->xMin > leaf.X) {
                    current->xMin = leaf.X;
                    return;
                }

                if (current->yMax < leaf.Y) {
                    current->yMax = leaf.Y;
                    return;
                }

                if (current->yMin > leaf.Y) {
                    current->yMin = leaf.Y;
                    return;
                }*/

                /*                                                  some prints for debug
                if (current->numNode >= 3){
                    int a = 0;
                    while (a < current->numNode){
                        printf("id - %d\n", ((leafNode*)current->node)[a].Hit.id);
                        a++;
                    }
                }
                 */
            }
        }else{          //it is not leaf node
            double newArea = DBL_MAX;
            int insertIndex = 0;
            for (int i = 0; i < current->numNode; i++){ //iterate over all nodes
                double nMaxX = ((Node*) current->node)[i].xMax, nMaxY = ((Node*) current->node)[i].yMax, nMinX = ((Node*) current->node)[i].xMin, nMinY = ((Node*) current->node)[i].yMin;

                if (((Node*) current->node)[i].xMin > leaf.X){  //compare the new point with every rectangle to know how much will increase the new area
                    nMinX = leaf.X;
                }

                if (((Node*) current->node)[i].yMin > leaf.Y){
                    nMinY = leaf.Y;
                }

                if (((Node*) current->node)[i].xMax < leaf.X){
                    nMaxX = leaf.X;
                }

                if (((Node*) current->node)[i].yMax < leaf.Y){
                    nMaxY = leaf.Y;
                }

                //calculate the are increased = new total area - old area
                double newAreaCreated = (nMaxY-nMinY)*(nMaxX-nMinX) - (((Node*) current->node)[i].xMax - ((Node*) current->node)[i].xMin)*(((Node*) current->node)[i].yMax - ((Node*) current->node)[i].yMin) ;
                if (newAreaCreated < newArea){  //insert in the minimum area created
                    newArea = newAreaCreated;
                    insertIndex = i;
                }else if (newAreaCreated == newArea){ //if equals insert in the smallest old area
                    if ((((Node*) current->node)[i].xMax - ((Node*) current->node)[i].xMin)*(((Node*) current->node)[i].yMax - ((Node*) current->node)[i].yMin) < (((Node*) current->node)[insertIndex].xMax - ((Node*) current->node)[insertIndex].xMin)*(((Node*) current->node)[insertIndex].yMax - ((Node*) current->node)[insertIndex].yMin)){
                        insertIndex = i;
                    }
                }
            }
            insert(leaf, &(((Node*) current->node)[insertIndex]));
            //back propagation

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
