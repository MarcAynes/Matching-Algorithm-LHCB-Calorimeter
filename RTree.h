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
    Node* parentNode;
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
                    current->numNode++;
                    ((leafNode*) (current->node))[0].Hit = leaf.Hit;
                    ((leafNode*) (current->node))[0].X = leaf.X;
                    ((leafNode*) (current->node))[0].Y = leaf.Y;
                    current->yMin = leaf.Y; //adding X and Y data
                    current->yMax = leaf.Y;
                    current->xMin = leaf.X;
                    current->xMax = leaf.X;
                }else {                     //does have info inserted
                    current->numNode++;
                    //copying leaf data to the final leafNode
                    ((leafNode *) (current->node))[(current->numNode) - 1].Hit = leaf.Hit;
                    ((leafNode *) (current->node))[(current->numNode) - 1].X = leaf.X;
                    ((leafNode *) (current->node))[(current->numNode) - 1].Y = leaf.Y;

                    if (current->xMax < leaf.X) {   //modifying X and Y data if needed
                        current->xMax = leaf.X;
                    }

                    if (current->xMin > leaf.X) {
                        current->xMin = leaf.X;
                    }

                    if (current->yMax < leaf.Y) {
                        current->yMax = leaf.Y;
                    }

                    if (current->yMin > leaf.Y) {
                        current->yMin = leaf.Y;
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
            }else{      // there is no space
                //adding leaf to the 6 temporally position we insert it to iterate over it.
                ((leafNode*) current->node)[current->numNode].Hit = leaf.Hit;
                ((leafNode*) current->node)[current->numNode].X = leaf.X;
                ((leafNode*) current->node)[current->numNode].Y = leaf.Y;
                int index1 = 0, index2 = 0;

                //find the 2 farthest points in the rectangle
                double dist = DBL_MAX;
                for (int i = 0; i <= current->numNode; i++) {     //O(N^2) cost at finding 2 furthest points, but for 6 points it's not critical
                    for (int j = i + 1; j <= current->numNode; j++) {
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

                //find the 2 elements closest to index1 both used to create the minimum rectangle
                double dist1 = DBL_MAX;
                double dist2 = DBL_MAX;
                int ind1 = 0, ind2 = 0;

                for (int i = 0; i <= current->numNode; i++){ //linear search for 2 elements closest to index1
                    if (i == index1 || i == index2){
                        continue;
                    }

                    double auxDist = euclideanDistance(((leafNode*)current->node)[i].X, ((leafNode*)current->node)[i].Y, ((leafNode*)current->node)[index1].X, ((leafNode*)current->node)[index1].Y);
                    if (auxDist < dist1){   //point i is nearest than previous one (1)
                        if (dist1 < dist2){ //point (1) is nearest than (2)
                            dist2 = dist1;
                            ind2 = ind1;
                        }
                        dist1 = auxDist;
                        ind1 = i;
                        continue;
                    }

                    if (auxDist < dist2){
                        dist2 = auxDist;
                        ind2 = i;
                    }
                }

                if (current->parentNode == nullptr){ // if the Node it's root
                    leafNode* auxNodes = (leafNode*) current->node;
                    current->node = nullptr;
                    current->node = (Node*) malloc(sizeof(Node)*(M+1)); //split, inserting 2 Nodes alloc all M+1 nodes because if we doa realloc we can lose the parent*
                    current->numNode = 2;
                    current->leaf = false;

                    //find the 2 other points which they will make the second rectangle
                    int ind3 = -1, ind4 = -1;
                    for (int i = 0; i < M + 1; i++){
                        if ((i == index1 || i == index2 || i == ind2 || i == ind1)) {
                            if (ind3 < 0) {
                                ind3 = i;
                            }else{
                                ind4 = i;
                                break;
                            }
                        }
                    }

                    //for each new rectangle insert m (3) leaf nodes
                    for (int j = 0; j < current->numNode; j++) {
                        //initializing each new Node
                        ((Node*) current->node)[j].yMin = DBL_MIN;
                        ((Node*) current->node)[j].xMin = DBL_MIN;
                        ((Node*) current->node)[j].xMax = DBL_MAX;
                        ((Node*) current->node)[j].yMax = DBL_MAX;
                        ((Node*) current->node)[j].node = (leafNode*) malloc(sizeof(leafNode)*(M + 1));
                        ((Node*) current->node)[j].numNode = 0;
                        ((Node*) current->node)[j].leaf = true;
                        ((Node*) current->node)[j].parentNode = current;

                        //adding working data of the leaf nodes
                        for (int i = 0; i < m; i++) {
                            int ind;
                            switch (j) {
                                case 0:
                                    switch (i) {
                                        case 0:
                                            ind = index1;
                                            break;
                                        case 1:
                                            ind = ind1;
                                            break;
                                        case 2:
                                            ind = ind2;
                                            break;
                                        default:
                                            ind = index1;
                                            break;
                                    }
                                    break;
                                default:
                                    switch (i) {
                                        case 0:
                                            ind = index2;
                                            break;
                                        case 1:
                                            ind = ind3;
                                            break;
                                        case 2:
                                            ind = ind4;
                                            break;
                                        default:
                                            ind = index2;
                                            break;
                                    }
                                    break;

                            }

                            //copying each leaf node to it's rectangle (Node)
                            ((Node*) current->node)[j].numNode++; //rectangle have +1 leaf node (data point)
                            ((leafNode*) ((Node*) current->node)[j].node)[((((Node*) current->node)[i].numNode) - 1)].Hit = auxNodes[ind].Hit; //adding hit data to it's rectangle (Node)
                            ((leafNode*) ((Node*) current->node)[j].node)[((((Node*) current->node)[i].numNode) - 1)].X = auxNodes[ind].X;
                            ((leafNode*) ((Node*) current->node)[j].node)[((((Node*) current->node)[i].numNode) - 1)].Y = auxNodes[ind].Y;

                            //modifying rectangle boundaries if needed
                            if (((Node*) current->node)[j].xMax < auxNodes[ind].X) {
                                ((Node*) current->node)[j].xMax = auxNodes[ind].X;
                            }

                            if (((Node*) current->node)[j].xMin > auxNodes[ind].X) {
                                ((Node*) current->node)[j].xMin = auxNodes[ind].X;
                            }

                            if (((Node*) current->node)[j].yMax < auxNodes[ind].Y) {
                                ((Node*) current->node)[j].yMax = auxNodes[ind].Y;
                            }

                            if (((Node*) current->node)[j].yMin > auxNodes[ind].Y) {
                                ((Node*) current->node)[j].yMin = auxNodes[ind].Y;
                            }
                        }
                    }
                }else{  //node is not root

                }
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
        root.node = (leafNode*) malloc(sizeof(leafNode)*(M + 1));   //alloc all memory, if we don't do this when we do a realloc on a Node we may lose the parent* node
        root.numNode = 0;
        root.xMax = -1;
        root.xMin = -1;
        root.yMax = -1;
        root.yMin = -1;
        root.parentNode = nullptr;
    }
};


#endif //MATCHINGALGORITHM_RTREE_H
