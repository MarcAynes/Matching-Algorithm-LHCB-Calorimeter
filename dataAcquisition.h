#ifndef MATCHINGALGORITHM_DATAACQUISITION_H
#define MATCHINGALGORITHM_DATAACQUISITION_H
#include <fstream>
#include <sstream>
#include <vector>
#include "hit.h"

using namespace std;

class dataAcquisition {
public:
    vector<vector<hit>> hitsData;

    /*
     * Reeding CSV file.
     * Storing hits data on vectors.
     */
    void readFile(string fileName){
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

                while(getline(str, word, ',')) {
                    if (word[0] == 'C' && element == 0) {
                        element++;

                    } else if(element == 0) {
                        i++;
                        break;
                    }

                    if(element == 1){
                        element++;
                        hitData.X = stod(word);
                    }

                    if(element == 2){
                        element++;
                        hitData.Y = stod(word);
                    }

                    if(element == 3){
                        element++;
                        hitData.Z = stod(word);
                    }

                    if(element == 4){
                        element++;
                        hitData.MeVs = stod(word);
                        data[i].push_back(hitData);
                    }
                }
            }
        }

        hitsData = data;
    }


};


#endif //MATCHINGALGORITHM_DATAACQUISITION_H
