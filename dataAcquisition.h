#ifndef MATCHINGALGORITHM_DATAACQUISITION_H
#define MATCHINGALGORITHM_DATAACQUISITION_H
#include <fstream>
#include <sstream>
#include <vector>
#include "hit.h"
#include "trace.h"

using namespace std;
/*
 * Reeding CSV file.
 * Storing hits data on vectors.
 */

vector<vector<hit>> readFile(string fileName, vector<vector<trace>>* traces){
    string line, word;
    fstream file (fileName, ios::in);
    vector<vector<hit>> data;
    int i = 0, k = 0;
    data.emplace_back();
    traces->emplace_back();

    if(file.is_open()){
        while(getline(file, line)){
            stringstream str(line);

            int element = 0;
            int j = 0;
            hit hitData;
            trace traceData;

            while(getline(str, word, ',')) {        //read line from csv
                if(element == 5){
                    element = 0;
                    hitData.MeVs = stod(word);
                    data[i].push_back(hitData);
                    j = 1;
                    break;
                }

                if(element == 4){
                    element++;
                    hitData.Z = stod(word);
                    continue;
                }

                if(element == 3){
                    element++;
                    hitData.Y = stod(word);
                    continue;
                }

                if(element == 2){
                    element++;
                    hitData.X = stod(word);
                    continue;
                }

                if (element == 1){
                    element++;
                    hitData.id = stoi(word);
                    continue;
                }

                if (word[0] == 'C' && element == 0) {
                    element++;
                } else if(element == 0) {
                    break;
                }
            }

            if (j == 1){
                continue;
            }

             do {        //read line from csv
                if (element == 6){
                    element++;
                    traceData.ty = stod(word);
                    (*traces)[i].push_back(traceData);
                    j = 1;
                    continue;
                }

                if (element == 5) {
                    element++;
                    traceData.tx = stod(word);
                    continue;
                }

                if (element == 4) {
                    element++;
                    traceData.Z = stod(word);
                    continue;
                }

                if (element == 3) {
                    element++;
                    traceData.Y = stod(word);
                    continue;
                }

                if (element == 2) {
                    element++;
                    traceData.X = stod(word);
                    continue;
                }

                if (element == 1) {
                    element++;
                    traceData.id = stoi(word);
                    continue;
                }

                if (word[0] == 'T' && element == 0) {
                    element++;
                    k = 0;
                } else if (element == 0 && k == 0) {
                    i++;
                    k++;
                    data.emplace_back();
                    traces->emplace_back();
                    break;
                }else if(k > 0){
                    break;
                }
            } while(getline(str, word, ','));

             if (j == 1){
                 continue;
             }
        }
    }

    return data;
}


#endif //MATCHINGALGORITHM_DATAACQUISITION_H
