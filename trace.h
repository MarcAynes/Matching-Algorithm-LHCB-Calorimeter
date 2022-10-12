//
// Created by Marc on 17/03/2022.
//

#ifndef MATCHINGALGORITHM_TRACE_H
#define MATCHINGALGORITHM_TRACE_H

#define Threshold 1000 //threshold in millimeters
#define MaxTraceHigh 12780 //Maximum Z where we will calculate the trace path

class trace {
public:
    double X;
    double Y;
    double Z;
    double tx;
    double ty;
    int id;
};


#endif //MATCHINGALGORITHM_TRACE_H
