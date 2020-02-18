//
// Created by ankur on 2/14/20.
//

#include "kmeans.h"
#include <stdlib.h>
#include <math.h>
#include <vector>

//#include <iostream>
//using namespace std;

int main() {
  int noOfItems;
  int k;

  noOfItems = 100;
  k = 10;

  int cluster[k];
  int oldCluster[k];
  int objects[noOfItems];
  int row[k];
  vector<vector<int>> groups;

  for (int i = 0 ; i < noOfItems ; ++i) {
    objects[i] = rand() % 50;
    if (i < k) {
      cluster[i] = objects[i]
    }
  }

  for (int i = 0 ; i < k ; ++i) {
    vector<int> newGroup;
    groups.push_back(newGroup);
  }

  int iter = 1;
  do {
    for (int i = 0 ; i < noOfItems ; ++i) {
      for (int j = 0 ; j < k ; ++j) {
        row[j] = abs
      }
    }
  }
}
