#ifndef MATCHINGALGORITHM_RTREE_H
#define MATCHINGALGORITHM_RTREE_H

#include "hit.h"
#include <vector>
#include <cfloat>
#include <bits/stdc++.h>

using namespace std;

#define m 3
#define M 5

class RTree {
private:
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

    Node root{};

    typedef struct Point{ //trace rectangle corners
        double Px;
        double Py;
    }Point;

    /*       ntop
     *    A ___ B   mPerpendicular
     *    /   /
 nLeft   /   /    nRight
   mm   /   /    mm
     * /___/    mPerpendicular
     * C   D
     *  ndown
     *
     *   Or
     *A ___ B
     * |   |
     * |   |
     * |   |
     * |___|
     *C     D
     */
    typedef struct traceRectangle{
        Point A;
        Point B;
        Point C;
        Point D;
        double mm;
        double mPerpendicular;
        double nTop;
        double nDown;
        double nRight;
        double nLeft;
        double minXLeft;
        double minXRight;
        double maxXLeft;
        double maxXRight;
        double minYLeft;
        double minYRight;
        double maxYLeft;
        double maxYRight;
        bool validM;    //if the first line is vertical m would be infinity
    }traceRectangle;

    traceRectangle traceRectangle1;

    /*
     * Function that prepares the state for the search
     * this function creates the trace rectangle corners
     */
    void prepareSearch(trace t){
        traceRectangle1.validM = false;

        //calculate the corners of the trace rectangle projection over rtree plane
        Point A; //initial point of the trace
        A.Px = t.X;
        A.Py = t.Y;

        Point B; //latest point of the trace
        B.Px = t.tx*(MaxTraceHigh - t.Z) + t.X;
        B.Py = t.ty*(MaxTraceHigh - t.Z) + t.Y;
        //Point B - Point A = trace line vector, also these 2 points are the beginning and end of the rectangle

        //Calculating 4 rectangle corners if it's tilted
        /*
         *    A ___ B
         *    /   /
         *   /   /
         *  /   /
         * /___/
         * C   D
         */
        //if the trace went through the detector perpendicular or almost perpendicular or there is no difference in X
        if ((B.Px > A.Px && (B.Px - A.Px) > Threshold) || (B.Px < A.Px && (A.Px - B.Px) > Threshold)){
            //Y = mX + n
            //m = (Y-Yo)/(X-Xo) => t.ty/t.tx
            //M and m names are in use by the min and max number of elements in a Rtree Node, we will use the name mm
            double mm = t.ty/t.tx;

            //n = Y - mX
            double n = A.Py - (mm * A.Px);

            //n + threshold = Y - m(X+threshold)
            //n - threshold = Y - m(X-threshold)
            double nLeft = A.Py - (mm * (A.Px + Threshold));    //left
            double nRight = A.Py - (mm * (A.Px - Threshold));   //right

            //Perpendicular lines m = -1/m
            double mPerpendicular = (-1)/mm;

            //perpendicular lines must cross through points previously calculated A and B
            double nTop = A.Py - (mPerpendicular * B.Px);
            double nDown = B.Py - (mPerpendicular * A.Px);

            //Calculating intersections between 4 lines X = (n2-n1)/(m1-m2)
            //Corner A of previous schema
            Point finalA;
            finalA.Px = (nTop-nLeft)/(mm - mPerpendicular);
            finalA.Py = mPerpendicular*finalA.Px + nTop;

            //Corner B of previous schema
            Point finalB;
            finalB.Px = (nTop-nRight)/(mm - mPerpendicular);
            finalB.Py = mPerpendicular*finalB.Px + nTop;

            //Corner C of previous schema
            Point finalC;
            finalC.Px = (nDown-nLeft)/(mm - mPerpendicular);
            finalC.Py = mPerpendicular*finalC.Px + nDown;

            //Corner D of previous schema
            Point finalD;
            finalD.Px = (nDown-nRight)/(mm - mPerpendicular);
            finalD.Py = mPerpendicular*finalD.Px + nDown;

            //copy values
            traceRectangle1.mm = mm;
            traceRectangle1.mPerpendicular = mPerpendicular;
            traceRectangle1.nDown = nDown;
            traceRectangle1.nTop = nTop;
            traceRectangle1.nRight = nRight;
            traceRectangle1.nLeft = nLeft;
            traceRectangle1.A = finalA;
            traceRectangle1.B = finalB;
            traceRectangle1.C = finalC;
            traceRectangle1.D = finalD;
            traceRectangle1.validM = true;

            if(traceRectangle1.A.Px < traceRectangle1.C.Px){
                traceRectangle1.minXLeft = traceRectangle1.A.Px;
                traceRectangle1.maxXLeft = traceRectangle1.C.Px;
            }else{
                traceRectangle1.minXLeft = traceRectangle1.C.Px;
                traceRectangle1.maxXLeft = traceRectangle1.A.Px;
            }

            if(traceRectangle1.B.Px < traceRectangle1.D.Px){
                traceRectangle1.minXRight = traceRectangle1.B.Px;
                traceRectangle1.maxXRight = traceRectangle1.D.Px;
            }else{
                traceRectangle1.minXRight = traceRectangle1.D.Px;
                traceRectangle1.maxXRight = traceRectangle1.B.Px;
            }

            if(traceRectangle1.B.Py < traceRectangle1.D.Py){
                traceRectangle1.minYRight = traceRectangle1.B.Py;
                traceRectangle1.maxYRight = traceRectangle1.D.Py;
            }else{
                traceRectangle1.minYRight = traceRectangle1.D.Py;
                traceRectangle1.maxYRight = traceRectangle1.B.Py;
            }

            if(traceRectangle1.A.Py < traceRectangle1.C.Py){
                traceRectangle1.minYLeft = traceRectangle1.A.Py;
                traceRectangle1.maxYLeft = traceRectangle1.C.Py;
            }else{
                traceRectangle1.minYLeft = traceRectangle1.C.Py;
                traceRectangle1.maxYLeft = traceRectangle1.A.Py;
            }
        }else{
            /*
             * if it's not tilted we calculate it in a different way
             *A ___ B
             * |   |
             * |   |
             * |   |
             * |___|
             *C     D
             */

            traceRectangle1.validM = false;  //vertical line has m = infinity on y = mx + n or almost

            if ((B.Py > A.Py && (B.Py - A.Py) < Threshold) || (B.Py < A.Py && (A.Py - B.Py) < Threshold)){             //the trace went through the detector almost perpendicular
                if(B.Py > A.Py){
                    //Corner A of previous schema
                    traceRectangle1.A.Px = B.Px - Threshold;
                    traceRectangle1.A.Py = B.Py + Threshold;

                    //Corner B of previous schema
                    traceRectangle1.B.Px = B.Px + Threshold;
                    traceRectangle1.B.Py = B.Py + Threshold;

                    //Corner C of previous schema
                    traceRectangle1.C.Px = A.Px - Threshold;
                    traceRectangle1.C.Py = A.Py - Threshold;

                    //Corner D of previous schema
                    traceRectangle1.D.Px = A.Px + Threshold;
                    traceRectangle1.D.Py = A.Py - Threshold;

                    traceRectangle1.mPerpendicular = 0;
                    traceRectangle1.nTop = traceRectangle1.A.Py;
                    traceRectangle1.nDown = traceRectangle1.C.Py;
                    traceRectangle1.nLeft = traceRectangle1.A.Px;
                    traceRectangle1.nRight = traceRectangle1.B.Px;
                }else{
                    //Corner A of previous schema
                    traceRectangle1.A.Px = A.Px - Threshold;
                    traceRectangle1.A.Py = A.Py + Threshold;

                    //Corner B of previous schema
                    traceRectangle1.B.Px = A.Px + Threshold;
                    traceRectangle1.B.Py = A.Py + Threshold;

                    //Corner C of previous schema
                    traceRectangle1.C.Px = B.Px - Threshold;
                    traceRectangle1.C.Py = B.Py - Threshold;

                    //Corner D of previous schema
                    traceRectangle1.D.Px = B.Px + Threshold;
                    traceRectangle1.D.Py = B.Py - Threshold;

                    traceRectangle1.mPerpendicular = 0;
                    traceRectangle1.nTop = traceRectangle1.A.Py;
                    traceRectangle1.nDown = traceRectangle1.C.Py;
                    traceRectangle1.nLeft = traceRectangle1.A.Px;
                    traceRectangle1.nRight = traceRectangle1.B.Px;
                }
            }else{                                      //the trace projection is vertical
                if(B.Py > A.Py) {
                    //Corner A of previous schema
                    traceRectangle1.A.Px = B.Px - Threshold;
                    traceRectangle1.A.Py = B.Py;

                    //Corner B of previous schema
                    traceRectangle1.B.Px = B.Px + Threshold;
                    traceRectangle1.B.Py = B.Py;

                    //Corner C of previous schema
                    traceRectangle1.C.Px = A.Px - Threshold;
                    traceRectangle1.C.Py = A.Py;

                    //Corner D of previous schema
                    traceRectangle1.D.Px = A.Px + Threshold;
                    traceRectangle1.D.Py = A.Py;

                    traceRectangle1.mPerpendicular = 0;
                    traceRectangle1.nTop = B.Py;
                    traceRectangle1.nDown = A.Py;
                    traceRectangle1.nLeft = traceRectangle1.A.Px;
                    traceRectangle1.nRight = traceRectangle1.B.Px;
                }else{
                    //Corner A of previous schema
                    traceRectangle1.A.Px = A.Px - Threshold;
                    traceRectangle1.A.Py = A.Py;

                    //Corner B of previous schema
                    traceRectangle1.B.Px = A.Px + Threshold;
                    traceRectangle1.B.Py = A.Py;

                    //Corner C of previous schema
                    traceRectangle1.C.Px = B.Px - Threshold;
                    traceRectangle1.C.Py = B.Py;

                    //Corner D of previous schema
                    traceRectangle1.D.Px = B.Px + Threshold;
                    traceRectangle1.D.Py = B.Py;

                    traceRectangle1.mPerpendicular = 0;
                    traceRectangle1.nTop = A.Py;
                    traceRectangle1.nDown = B.Py;
                    traceRectangle1.nLeft = traceRectangle1.A.Px;
                    traceRectangle1.nRight = traceRectangle1.B.Px;
                }
            }

            if(traceRectangle1.B.Py < traceRectangle1.D.Py){
                traceRectangle1.minYRight = traceRectangle1.B.Py;
                traceRectangle1.maxYRight = traceRectangle1.D.Py;
            }else{
                traceRectangle1.minYRight = traceRectangle1.D.Py;
                traceRectangle1.maxYRight = traceRectangle1.B.Py;
            }

            if(traceRectangle1.A.Py < traceRectangle1.C.Py){
                traceRectangle1.minYLeft = traceRectangle1.A.Py;
                traceRectangle1.maxYLeft = traceRectangle1.C.Py;
            }else{
                traceRectangle1.minYLeft = traceRectangle1.C.Py;
                traceRectangle1.maxYLeft = traceRectangle1.A.Py;
            }

            if(traceRectangle1.B.Px < traceRectangle1.A.Px){
                traceRectangle1.minXRight = traceRectangle1.A.Px;
                traceRectangle1.minXRight = traceRectangle1.A.Px;
                traceRectangle1.maxXLeft = traceRectangle1.B.Px;
                traceRectangle1.minXLeft = traceRectangle1.B.Px;
            }else{
                traceRectangle1.minXRight = traceRectangle1.B.Px;
                traceRectangle1.maxXRight = traceRectangle1.B.Px;
                traceRectangle1.maxXLeft = traceRectangle1.A.Px;
                traceRectangle1.minXLeft = traceRectangle1.A.Px;
            }
        }
    }

    /*
     * Function that given a trace (mathematical vector) + threshold will create a rectangle
     * Given this rectangle will calculate if the leaf Node (Hit) is inside
     */
    bool hitIsInside(trace t, leafNode leaf){
        if (traceRectangle1.validM) {
            double n = leaf.Y - (traceRectangle1.mm * leaf.X);

            // y = ((n1*m2)-(n2*m1))/(m2-m1) <- used to know the Y of the trace path which must be between
            // 1 elements = leafNode    2 elements = trace
            double y1 = ((n * traceRectangle1.mPerpendicular) - (traceRectangle1.nDown * traceRectangle1.mm)) / (traceRectangle1.mPerpendicular - traceRectangle1.mm);
            double y2 = ((n * traceRectangle1.mPerpendicular) - (traceRectangle1.nTop * traceRectangle1.mm)) / (traceRectangle1.mPerpendicular - traceRectangle1.mm);

            if ((traceRectangle1.nLeft < traceRectangle1.nRight && n <= traceRectangle1.nRight && n >= traceRectangle1.nLeft || traceRectangle1.nLeft > traceRectangle1.nRight && n >= traceRectangle1.nRight && n <= traceRectangle1.nLeft) && leaf.Y >= y1 && leaf.Y <= y2) {
                return true;
            }
        }else{
            if (traceRectangle1.minXLeft <= leaf.X && traceRectangle1.maxXRight >= leaf.X && traceRectangle1.minYRight <= leaf.Y && traceRectangle1.maxYRight >= leaf.Y){
                return true;
            }
        }
        return false;
    }

    /*
     * Function that given a trace (mathematical vector) + threshold will create a rectangle
     * Given this rectangle will calculate if the Node is partially or totally inside
     */
    bool isInside(trace t, Node node){
        /* the rectangle (Node) have this shape
         *  ___
         * |   |
         * |___|
         *
         * 2 of the lines de m values on y = mx + n is 0
         * the other 2 m are infinity
         */
        double y;
        if (traceRectangle1.validM){    // the trace does not have a vertical path
            //y = mx + n -> using x as de vertical left side of the rectangle (Node)
            y = (traceRectangle1.mm * node.xMin) + traceRectangle1.nLeft;

            if (y <= node.yMax && node.yMin >= y && y <= traceRectangle1.maxYLeft && y >= traceRectangle1.minYLeft) {     //if Y is "touching" both rectangles
                return true;
            }

            //y = mx + n -> using x as de vertical right side of the rectangle (Node)
            y = (traceRectangle1.mm * node.xMax) + traceRectangle1.nLeft;

            if (y <= node.yMax && node.yMin >= y && y <= traceRectangle1.maxYLeft && y >= traceRectangle1.minYLeft) {     //if Y is "touching" both rectangles
                return true;
            }


            //y = mx + n -> using x as de vertical left side of the rectangle (Node)
            y = (traceRectangle1.mm * node.xMin) + traceRectangle1.nRight;

            if (y <= node.yMax && node.yMin >= y && y <= traceRectangle1.maxYRight && y >= traceRectangle1.minYRight) {     //if Y is "touching" both rectangles
                return true;
            }

            //y = mx + n -> using x as de vertical right side of the rectangle (Node)
            y = (traceRectangle1.mm * node.xMax) + traceRectangle1.nRight;

            if (y <= node.yMax && node.yMin >= y && y <= traceRectangle1.maxYRight && y >= traceRectangle1.minYRight) {     //if Y is "touching" both rectangles
                return true;
            }

            //check if there is an intersection between trace rectangle left or right and Node rectangle up or down sides
            double x;
            if(node.yMin >= traceRectangle1.minYRight && node.yMin <= traceRectangle1.maxYRight) {
                //node rectangle down line
                //(y-n)/m = x
                x = ((((double) node.yMin) - ((double) traceRectangle1.nRight)) / (double) traceRectangle1.mm);

                if (x <= node.xMax && node.xMin <= x){  //trace rectangle, right line
                    return true;
                }

                x = ((((double) node.yMin) - ((double) traceRectangle1.nLeft)) / (double) traceRectangle1.mm);

                if (x <= node.xMax && node.xMin <= x){ //trace rectangle, left line
                    return true;
                }
            }

            if(node.yMax >= traceRectangle1.minYRight && node.yMax <= traceRectangle1.maxYRight) {
                //node rectangle top line
                //(y-n)/m = x
                x = ((((double) node.yMax) - ((double) traceRectangle1.nRight)) / (double) traceRectangle1.mm);

                if (x <= node.xMax && node.xMin <= x){  //trace rectangle, right line
                    return true;
                }

                x = ((((double) node.yMax) - ((double) traceRectangle1.nLeft)) / (double) traceRectangle1.mm);

                if (x <= node.xMax && node.xMin <= x){ //trace rectangle, left line
                    return true;
                }
            }

            //if trace is inside a node rectangle, we will compare if one of the rectangle vertex is inside node
            if (traceRectangle1.A.Py <= node.yMax && traceRectangle1.A.Py >= node.yMin && traceRectangle1.A.Px <= node.xMax && traceRectangle1.A.Px >= node.xMin){
                return true;
            }

            //if node rectangle is inside trace rectangle, making a line through one corner of node rectangle with m = trace m
            // if trace nleft < n corner < trace nright => between corner is between 2 sides of trace rectangle, then is needed to know if the Y of the corner is between ymax and ymin of the trace
            //y-mx = n
            double n = node.yMin - (traceRectangle1.mm * node.xMin);

            // y = ((n1*m2)-(n2*m1))/(m2-m1) <- used to know the Y of the trace path which must be between
            // 1 elements = Node    2 elements = trace
            double y1 = ((n * traceRectangle1.mPerpendicular)-(traceRectangle1.nDown * traceRectangle1.mm))/(traceRectangle1.mPerpendicular - traceRectangle1.mm);
            double y2 = ((n * traceRectangle1.mPerpendicular)-(traceRectangle1.nTop * traceRectangle1.mm))/(traceRectangle1.mPerpendicular - traceRectangle1.mm);

            if ((traceRectangle1.nLeft < traceRectangle1.nRight && n <= traceRectangle1.nRight && n >= traceRectangle1.nLeft || traceRectangle1.nLeft > traceRectangle1.nRight && n >= traceRectangle1.nRight && n <= traceRectangle1.nLeft) && node.yMin >= y1 && node.yMin <= y2){
                return true;
            }

            return false;
        }else{ //trace have a vertical rectangle
            //rectangles are intersecting
            if (traceRectangle1.minXLeft <= node.xMax && traceRectangle1.maxXRight >= node.xMin && traceRectangle1.minYRight <= node.yMax && traceRectangle1.maxYRight >= node.yMin){
                return true;
            }
        }
        return false;
    }

    double euclideanDistance(double x1, double y1, double x2, double y2){
        // Calculating distance
        return sqrt(pow(x2 - x1, 2) +
                    pow(y2 - y1, 2));
    }

    /*
     * Returns: Area of the rectangle needed to insert 2 little rectangle minus both areas of the little rectangles
     */
    double extraArea(Node littleRectangleA, Node littleRectangleB){
        double xMin = littleRectangleA.xMin < littleRectangleB.xMin ? littleRectangleA.xMin : littleRectangleB.xMin;
        double yMin = littleRectangleA.yMin < littleRectangleB.yMin ? littleRectangleA.yMin : littleRectangleB.yMin;
        double xMax = littleRectangleA.xMax > littleRectangleB.xMax ? littleRectangleA.xMax : littleRectangleB.xMax;
        double yMax = littleRectangleA.yMax > littleRectangleB.yMax ? littleRectangleA.yMax : littleRectangleB.yMax;

        double totalArea = (xMax - xMin) * (yMax - yMin);
        double areaA = (littleRectangleA.xMax - littleRectangleA.xMin) * (littleRectangleA.yMax - littleRectangleA.yMin);
        double areaB = (littleRectangleB.xMax - littleRectangleB.xMin) * (littleRectangleB.yMax - littleRectangleB.yMin);

        return totalArea - areaA - areaB;
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
                current->numNode++;
                int index1 = 0, index2 = 0;

                //find the 2 farthest points in the rectangle
                double dist = DBL_MIN;
                for (int i = 0; i < current->numNode; i++) {     //O(N^2) cost at finding 2 furthest points, but for 6 points it's not critical
                    for (int j = i + 1; j < current->numNode; j++) {
                        double auxDist = euclideanDistance(((leafNode *) current->node)[i].X,
                                                           ((leafNode *) current->node)[i].Y,
                                                           ((leafNode *) current->node)[j].X,
                                                           ((leafNode *) current->node)[j].Y);
                        if (auxDist > dist) {
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

                for (int i = 0; i < current->numNode; i++){ //linear search for 2 elements closest to index1
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

                //find the 2 other points which they will make the second rectangle
                int ind3 = -1, ind4 = -1;
                for (int i = 0; i < M + 1; i++){
                    if ((i != index1 && i != index2 && i != ind2 && i != ind1)) {
                        if (ind3 < 0) {
                            ind3 = i;
                        }else{
                            ind4 = i;
                            break;
                        }
                    }
                }

                //It is not important if we insert in a full (5 nodes) parent a new node
                //due to when we exit this function, because it's recursive, we will be in the case
                //where we need to rotate nodes (split)
                leafNode* auxNodes = (leafNode*) current->node;
                current->node = nullptr;

                if (current->parentNode == nullptr){ // if the Node it's root
                    current->node = (Node*) malloc(sizeof(Node)*(M+1));
                    current->numNode = 2;
                    current->leaf = false;

                    //for each new rectangle insert m (3) leaf nodes
                    for (int j = 0; j < 2; j++){
                        //initializing each new Node
                        ((Node*) current->node)[j].yMin = (double) 999999999999999999;
                        ((Node*) current->node)[j].xMin = (double) 999999999999999999;  //DBL_MAX and DBL_MIN sometimes not working properly
                        ((Node*) current->node)[j].xMax = (double) -999999999999999999;
                        ((Node*) current->node)[j].yMax = (double) -999999999999999999;
                        ((Node*) current->node)[j].node = (leafNode*) malloc(sizeof(leafNode)*(M + 1));
                        ((Node*) current->node)[j].numNode = 0;
                        ((Node*) current->node)[j].leaf = true;
                        ((Node*) current->node)[j].parentNode = current;

                        for (int i = 0; i < m; i++) {
                            int ind = j == 0 ? i == 0 ? index1 : i == 1 ? ind1 : ind2     :     i == 0 ? index2 : i == 1 ? ind3 : ind4 ;

                            //copying each leaf node to it's rectangle (Node)
                            ((Node*) current->node)[j].numNode++; //rectangle have +1 leaf node (data point)
                            ((leafNode*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].Hit = auxNodes[ind].Hit; //adding hit data to it's rectangle (Node)
                            ((leafNode*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].X = auxNodes[ind].X;
                            ((leafNode*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].Y = auxNodes[ind].Y;

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

                    //modify current (root) rectangle size. Must include both of them
                    current->xMin = ((Node*) current->node)[0].xMin < ((Node*) current->node)[1].xMin ? ((Node*) current->node)[0].xMin : ((Node*) current->node)[1].xMin;
                    current->yMin = ((Node*) current->node)[0].yMin < ((Node*) current->node)[1].yMin ? ((Node*) current->node)[0].yMin : ((Node*) current->node)[1].yMin;
                    current->xMax = ((Node*) current->node)[0].xMax > ((Node*) current->node)[1].xMax ? ((Node*) current->node)[0].xMax : ((Node*) current->node)[1].xMax;
                    current->yMax = ((Node*) current->node)[0].yMax > ((Node*) current->node)[1].yMax ? ((Node*) current->node)[0].yMax : ((Node*) current->node)[1].yMax;

                }else{  //node is not root
                    //current->node = (leafNode *) malloc(sizeof(leafNode)*(M+1)); //split, inserting 2 Nodes alloc all M+1 nodes because if we do a realloc we can lose the parent*
                    current->parentNode->numNode++;         // this is important if numNodes > M we must split when we exit the function
                    current->numNode = 0;
                    //we have index to the current position and the new rectangle in the last parnet->node position created
                    int currentIndexInParent = (int) (current - ((Node*) current->parentNode->node));
                    //now we need to divide the M+1 leafs into those 2 nodes. m elements each one.

                    for (int x = 0; x < 2; x++){
                        //if we are root node we will only have 2 elements, so we need to equal x else we may have the index in different positions
                        int j = x == 0 ? currentIndexInParent : (current->parentNode->numNode - 1);
                        //initializing each new Node
                        ((Node*) current->parentNode->node)[j].yMin = (double) 999999999999999999;
                        ((Node*) current->parentNode->node)[j].xMin = (double) 999999999999999999;
                        ((Node*) current->parentNode->node)[j].xMax = (double) -999999999999999999;
                        ((Node*) current->parentNode->node)[j].yMax = (double) -999999999999999999;
                        ((Node*) current->parentNode->node)[j].node = (leafNode*) malloc(sizeof(leafNode)*(M + 1));
                        ((Node*) current->parentNode->node)[j].numNode = 0;
                        ((Node*) current->parentNode->node)[j].leaf = true;
                        ((Node*) current->parentNode->node)[j].parentNode = current->parentNode;

                        for (int i = 0; i < m; i++) {
                            //different index depends on Node where are we working (which of 2 rectangles we are inserting data)
                            int ind = j == currentIndexInParent ? i == 0 ? index1 : i == 1 ? ind1 : ind2     :     i == 0 ? index2 : i == 1 ? ind3 : ind4 ;

                            //copying each leaf node to it's rectangle (Node)
                            ((Node*) current->parentNode->node)[j].numNode++; //rectangle have +1 leaf node (data point)
                            ((leafNode*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].Hit = auxNodes[ind].Hit; //adding hit data to it's rectangle (Node)
                            ((leafNode*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].X = auxNodes[ind].X;
                            ((leafNode*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].Y = auxNodes[ind].Y;

                            //modifying rectangle boundaries if needed
                            if (((Node*) current->parentNode->node)[j].xMax < auxNodes[ind].X) {
                                ((Node*) current->parentNode->node)[j].xMax = auxNodes[ind].X;
                            }

                            if (((Node*) current->parentNode->node)[j].xMin > auxNodes[ind].X) {
                                ((Node*) current->parentNode->node)[j].xMin = auxNodes[ind].X;
                            }

                            if (((Node*) current->parentNode->node)[j].yMax < auxNodes[ind].Y) {
                                ((Node*) current->parentNode->node)[j].yMax = auxNodes[ind].Y;
                            }

                            if (((Node*) current->parentNode->node)[j].yMin > auxNodes[ind].Y) {
                                ((Node*) current->parentNode->node)[j].yMin = auxNodes[ind].Y;
                            }
                        }
                    }
                    //parent size will change when we exit this recursive call. Every current node should modify himself
                }
                free(auxNodes); //free copied data
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
            insert(leaf, &(((Node*) current->node)[insertIndex]));  //recursive call
            //back propagation
            //modify X and Y limits if leaf node it's outside the actual rectangle
            current->xMin = current->xMin > leaf.X ? leaf.X : current->xMin;
            current->yMin = current->yMin > leaf.Y ? leaf.Y : current->yMin;
            current->xMax = current->xMax < leaf.X ? leaf.X : current->xMax;
            current->yMax = current->yMax < leaf.Y ? leaf.Y : current->yMax;

            if (current->numNode > M){ //we have M+1 Nodes, we need to split
                int currentIndexInParent;

                double area = DBL_MIN;
                int index1 = 0, index2 = 0;
                for (int i = 0; i < current->numNode; i++) {     //O(N^2) cost at finding 2 furthest rectangles, but for M + 1 = 6 rectangles it's not critical
                    for (int j = i + 1; j < current->numNode; j++) {
                        double auxArea = extraArea(((Node*) current->node)[i], ((Node*) current->node)[j]);
                        if (auxArea > area) {
                            index1 = i;
                            index2 = j;
                            area = auxArea;
                        }
                    }
                }

                //find the 2 rectangles closest to index1 both used to create the minimum rectangle possible
                double area1 = DBL_MAX;
                double area2 = DBL_MAX;
                int ind1 = 0, ind2 = 0;

                for (int i = 0; i < current->numNode; i++){ //linear search for 2 elements closest to index1
                    if (i == index1 || i == index2){
                        continue;
                    }

                    double auxArea = extraArea(((Node *)current->node)[i], ((Node*)current->node)[index1]);
                    if (auxArea < area1){   //rectangle i is nearest than previous one (1)
                        if (area1 < area2){ //rectangle (1) is nearest than (2)
                            area2 = area1;
                            ind2 = ind1;
                        }
                        area1 = auxArea;
                        ind1 = i;
                        continue;
                    }

                    if (auxArea < area2){
                        area2 = auxArea;
                        ind2 = i;
                    }
                }

                //find the 2 other points which they will make the second rectangle
                int ind3 = -1, ind4 = -1;
                for (int i = 0; i < M + 1; i++){
                    if ((i != index1 && i != index2 && i != ind2 && i != ind1)) {
                        if (ind3 < 0) {
                            ind3 = i;
                        }else{
                            ind4 = i;
                            break;
                        }
                    }
                }
                Node* auxNodes = (Node*) current->node;
                current->node = nullptr;

                if (current->parentNode == nullptr){ // root node splitting root node in 2, creating new nodes
                    current->node = (Node*) malloc(sizeof(Node)*(M+1));
                    current->numNode = 2;
                    current->leaf = false;

                    for (int j = 0; j < 2; j++){
                        //initializing each new Node
                        ((Node*) current->node)[j].yMin = (double) 999999999999999999;
                        ((Node*) current->node)[j].xMin = (double) 999999999999999999;
                        ((Node*) current->node)[j].xMax = (double) -999999999999999999;
                        ((Node*) current->node)[j].yMax = (double) -999999999999999999;
                        ((Node*) current->node)[j].node = (Node *) malloc(sizeof(Node)*(M + 1));
                        ((Node*) current->node)[j].numNode = 0;
                        ((Node*) current->node)[j].leaf = false;
                        ((Node*) current->node)[j].parentNode = current;

                        for (int i = 0; i < m; i++) {
                            int ind = j == 0 ? i == 0 ? index1 : i == 1 ? ind1 : ind2     :     i == 0 ? index2 : i == 1 ? ind3 : ind4 ;

                            //copying each child Node to its new parent Node (Node)
                            ((Node*) current->node)[j].numNode++; //rectangle have +1 leaf node (data point)
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].node = auxNodes[ind].node; //adding Node data to its parent rectangle (Node)
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].numNode = auxNodes[ind].numNode;
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].xMax = auxNodes[ind].xMax;
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].yMin = auxNodes[ind].yMin;
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].xMin = auxNodes[ind].xMin;
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].yMax = auxNodes[ind].yMax;
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].leaf = auxNodes[ind].leaf;
                            ((Node*) ((Node*) current->node)[j].node)[((((Node*) current->node)[j].numNode) - 1)].parentNode = &((Node*)(current->node))[j];

                            //modifying rectangle boundaries if needed
                            if (((Node*) current->node)[j].xMax < auxNodes[ind].xMax) {
                                ((Node*) current->node)[j].xMax = auxNodes[ind].xMax;
                            }

                            if (((Node*) current->node)[j].xMin > auxNodes[ind].xMin) {
                                ((Node*) current->node)[j].xMin = auxNodes[ind].xMin;
                            }

                            if (((Node*) current->node)[j].yMax < auxNodes[ind].yMax) {
                                ((Node*) current->node)[j].yMax = auxNodes[ind].yMax;
                            }

                            if (((Node*) current->node)[j].yMin > auxNodes[ind].yMin) {
                                ((Node*) current->node)[j].yMin = auxNodes[ind].yMin;
                            }
                        }
                    }
                    //modify current (root) rectangle size. Must include both new nodes
                    current->xMin = ((Node*) current->node)[0].xMin < ((Node*) current->node)[1].xMin ? ((Node*) current->node)[0].xMin : ((Node*) current->node)[1].xMin;
                    current->yMin = ((Node*) current->node)[0].yMin < ((Node*) current->node)[1].yMin ? ((Node*) current->node)[0].yMin : ((Node*) current->node)[1].yMin;
                    current->xMax = ((Node*) current->node)[0].xMax > ((Node*) current->node)[1].xMax ? ((Node*) current->node)[0].xMax : ((Node*) current->node)[1].xMax;
                    current->yMax = ((Node*) current->node)[0].yMax > ((Node*) current->node)[1].yMax ? ((Node*) current->node)[0].yMax : ((Node*) current->node)[1].yMax;

                }else{ //not root node
                    //split, inserting 2 Nodes alloc all M+1 nodes
                    current->parentNode->numNode++;         // this is important if numNodes > M we must split when we exit the function (back propagation)
                    current->numNode = 0;
                    //we have an index to the current position and the new rectangle in the last parnet->node position created
                    int currentIndexInParent = (int) (current - ((Node*) current->parentNode->node));
                    //now we need to divide the M+1 leafs into those 2 nodes. m elements each one.

                    for (int x = 0; x < 2; x++){
                        //in the first iteration we will be inserting m nodes to the current node, in the second one we will insert it in the new node.
                        int j = x == 0 ? currentIndexInParent : (current->parentNode->numNode - 1);
                        //initializing each new Node
                        ((Node*) current->parentNode->node)[j].yMin = (double) 999999999999999999;
                        ((Node*) current->parentNode->node)[j].xMin = (double) 999999999999999999;
                        ((Node*) current->parentNode->node)[j].xMax = (double) -999999999999999999;
                        ((Node*) current->parentNode->node)[j].yMax = (double) -999999999999999999;
                        ((Node*) current->parentNode->node)[j].node = (Node*) malloc(sizeof(Node)*(M + 1));
                        ((Node*) current->parentNode->node)[j].numNode = 0;
                        ((Node*) current->parentNode->node)[j].leaf = false;
                        ((Node*) current->parentNode->node)[j].parentNode = current->parentNode;

                        for (int i = 0; i < m; i++) {
                            //different index depends on Node where are we working (which of 2 rectangles we are inserting data)
                            int ind = j == currentIndexInParent ? i == 0 ? index1 : i == 1 ? ind1 : ind2     :     i == 0 ? index2 : i == 1 ? ind3 : ind4 ;

                            //copying each leaf node to its rectangle (Node)
                            ((Node*) current->parentNode->node)[j].numNode++; //rectangle have +1 leaf node (data point)
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].node = auxNodes[ind].node; //adding Node data to its parent rectangle (Node)
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].numNode = auxNodes[ind].numNode;
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].xMax = auxNodes[ind].xMax;
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].yMin = auxNodes[ind].yMin;
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].xMin = auxNodes[ind].xMin;
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].yMax = auxNodes[ind].yMax;
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].leaf = auxNodes[ind].node;
                            ((Node*) ((Node*) current->parentNode->node)[j].node)[((((Node*) current->parentNode->node)[j].numNode) - 1)].parentNode = &(((Node*) current->parentNode->node)[j]);

                            //modifying current and new rectangle boundaries
                            if (((Node*) current->parentNode->node)[j].xMax < auxNodes[ind].xMax) {
                                ((Node*) current->parentNode->node)[j].xMax = auxNodes[ind].xMax;
                            }

                            if (((Node*) current->parentNode->node)[j].xMin > auxNodes[ind].xMin) {
                                ((Node*) current->parentNode->node)[j].xMin = auxNodes[ind].xMin;
                            }

                            if (((Node*) current->parentNode->node)[j].yMax < auxNodes[ind].yMax) {
                                ((Node*) current->parentNode->node)[j].yMax = auxNodes[ind].yMax;
                            }

                            if (((Node*) current->parentNode->node)[j].yMin > auxNodes[ind].yMin) {
                                ((Node*) current->parentNode->node)[j].yMin = auxNodes[ind].yMin;
                            }
                        }
                    }
                }
                free(auxNodes);
            }
        }
    }

    /*
     * Recursive function where we will search for a given trace.
     */
    void searchTrace(trace traceIn, Node *current, vector<hit> *data){
        if (current->leaf){
            for (int i = 0; i < current->numNode; i++){
                if(hitIsInside(traceIn, ((leafNode*)(current->node))[i])){
                    data->push_back(((leafNode*)(current->node))[i].Hit);
                }
            }
        }else{
            for (int i = 0; i < current->numNode; i++){
                if (isInside(traceIn, ((Node*)(current->node))[i])){
                    searchTrace(traceIn, &((Node*)(current->node))[i], data);
                }
            }
        }
    }

public:
    /*
     * Insert Function
     * add ALL data to the data structure,
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

    /*
     * Search ALL possible hits inside a trace path + some threshold
     * Using recursive function we iterate over all the RTree be aware of the stack used and the stack available
     */
    vector<vector<vector<hit>>> searchWithTrace(vector<vector<trace>> traces){
        vector<vector<vector<hit>>> possibleHits;
        possibleHits.emplace_back();
        for (int i = 0; i < traces.size(); i++) {
            possibleHits[i].emplace_back();
            for (int j = 0; j < traces.size(); j++){
                possibleHits[i][j].emplace_back();
                vector<hit> data;
                data.emplace_back();
                prepareSearch(traces[i][j]);
                searchTrace(traces[i][j], &root, &data);
                possibleHits[i][j] = data;
            }
        }
        return possibleHits;
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
