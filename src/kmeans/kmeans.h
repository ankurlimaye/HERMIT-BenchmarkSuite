// kmeans.h
// Ethan Brodsky
// October 2011

//
// Created by ankur on 2/14/20.
//

#ifndef KMEANS_H_
#define KMEANS_H_

void kmeans(
  int dim,                        // dimension of data
  double *X,                      // pointer to data
  int n,                          // number of elements
  int k,                          // number of clusters
  double *cluster_centroid,       // initial cluster centroids
  int *cluster_assignment_final   // output
);

#endif //KMEANS_H_
