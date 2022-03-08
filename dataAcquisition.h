#ifndef MATCHINGALGORITHM_DATAACQUISITION_H
#define MATCHINGALGORITHM_DATAACQUISITION_H
#include <fstream>
#include <sstream>
#include <vector>
#include "hit.h"

using namespace std;
/*
 * Reeding CSV file.
 * Storing hits data on vectors.
 */

vector<vector<hit>> readFile(string fileName){
    string line, word;
    fstream file (fileName, ios::in);
    vector<vector<hit>> data;
    int i = 0;

    if(file.is_open()){
        while(getline(file, line)){
            stringstream str(line);

            int element = 0;
            int j = 0;
            hit hitData;

            data.emplace_back();

            while(getline(str, word, ',')) {        //read line from csv
                if(element == 5){
                    element++;
                    hitData.MeVs = stod(word);
                    data[i].push_back(hitData);
                    continue;
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
                    i++;
                    break;
                }
            }
        }
    }

    return data;
}


#endif //MATCHINGALGORITHM_DATAACQUISITION_H
