// kmeans.c
// Ethan Brodsky
// October 2011

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#define sqr(x) ((x)*(x))
#define MAX_CLUSTERS 10
#define MAX_ITERATIONS 100
#define BIG_double (INFINITY)

void fail(char *str) {
  printf("%s", str);
  exit(-1);
}

double calc_distance(int dim, double *p1, double *p2) {
  double distance_sq_sum = 0;

  for (int ii = 0 ; ii < dim ; ii++) {
    distance_sq_sum += sqr(p1[ii] - p2[ii]);
  }

  return distance_sq_sum;
}

void calc_all_distances(int dim, int n, int k, double *X, double *centroid, double *distance_output) {
  for (int ii = 0 ; ii < n ; ii++) { // for each point
    for (int jj = 0 ; jj < k ; jj++) { // for each cluster
      // calculate distance between point and cluster centroid
      distance_output[ii * k + jj] = calc_distance(dim, &X[ii * dim], &centroid[jj * dim]);
    }
  }
}

double calc_total_distance(int dim, int n, int k, double *X, double *centroids, int *cluster_assignment_index) {
// NOTE: a point with cluster assignment -1 is ignored

  double tot_D = 0;

  // for every point
  for (int ii = 0 ; ii < n ; ii++) {
    // which cluster is it in?
    int active_cluster = cluster_assignment_index[ii];

    // sum distance
    if (active_cluster != -1) {
      tot_D += calc_distance(dim, &X[ii * dim], &centroids[active_cluster * dim]);
    }
  }

  return tot_D;
}

void choose_all_clusters_from_distances(int dim, int n, int k, double *distance_array, int *cluster_assignment_index) {
  // for each point
  for (int ii = 0 ; ii < n ; ii++) {
    int best_index = -1;
    double closest_distance = BIG_double;

    // for each cluster
    for (int jj = 0 ; jj < k ; jj++) {
      // distance between point and cluster centroid

      double cur_distance = distance_array[ii * k + jj];
      if (cur_distance < closest_distance) {
        best_index = jj;
        closest_distance = cur_distance;
      }
    }

    // record in array
    cluster_assignment_index[ii] = best_index;
  }
}

void calc_cluster_centroids(int dim, int n, int k, double *X, int *cluster_assignment_index, double *new_cluster_centroid) {
  int cluster_member_count[MAX_CLUSTERS];

  // initialize cluster centroid coordinate sums to zero
  for (int ii = 0 ; ii < k ; ii++) {
    cluster_member_count[ii] = 0;

    for (int jj = 0 ; jj < dim ; jj++) {
      new_cluster_centroid[ii * dim + jj] = 0;
    }
  }

  // sum all points
  // for every point
  for (int ii = 0 ; ii < n ; ii++) {
    // which cluster is it in?
    int active_cluster = cluster_assignment_index[ii];

    // update count of members in that cluster
    cluster_member_count[active_cluster]++;

    // sum point coordinates for finding centroid
    for (int jj = 0 ; jj < dim ; jj++) {
      new_cluster_centroid[active_cluster * dim + jj] += X[ii * dim + jj];
    }
  }

  // now divide each coordinate sum by number of members to find mean/centroid
  // for each cluster
  for (int ii = 0 ; ii < k ; ii++) {
    if (cluster_member_count[ii] == 0) {
      printf("WARNING: Empty cluster %d! \n", ii);
    }

    // for each dimension
    for (int jj = 0 ; jj < dim ; jj++) {
      new_cluster_centroid[ii * dim + jj] /= cluster_member_count[ii];  /// XXXX will divide by zero here for any empty clusters!
    }
  }
}

void get_cluster_member_count(int n, int k, int *cluster_assignment_index, int *cluster_member_count) {
  // initialize cluster member counts
  for (int ii = 0 ; ii < k ; ii++) {
    cluster_member_count[ii] = 0;
  }

  // count members of each cluster
  for (int ii = 0 ; ii < n ; ii++) {
    cluster_member_count[cluster_assignment_index[ii]]++;
  }
}

void update_delta_score_table(int dim, int n, int k, double *X, int *cluster_assignment_cur, double *cluster_centroid, int *cluster_member_count, double *point_move_score_table, int cc) {
  // for every point (both in and not in the cluster)
  for (int ii = 0 ; ii < n ; ii++) {
    double dist_sum = 0;
    for (int kk = 0 ; kk < dim ; kk++) {
      double axis_dist = X[ii * dim + kk] - cluster_centroid[cc * dim + kk];
      dist_sum += sqr(axis_dist);
    }

    double mult = ((double) cluster_member_count[cc] / (cluster_member_count[cc] + ((cluster_assignment_cur[ii] == cc) ? -1 : +1)));

    point_move_score_table[ii * dim + cc] = dist_sum * mult;
  }
}

void perform_move(int dim, int n, int k, double *X, int *cluster_assignment, double *cluster_centroid, int *cluster_member_count, int move_point, int move_target_cluster) {
  int cluster_old = cluster_assignment[move_point];
  int cluster_new = move_target_cluster;

  // update cluster assignment array
  cluster_assignment[move_point] = cluster_new;

  // update cluster count array
  cluster_member_count[cluster_old]--;
  cluster_member_count[cluster_new]++;

  if (cluster_member_count[cluster_old] <= 1)
    printf("WARNING: Can't handle single-member clusters! \n");

  // update centroid array
  for (int ii = 0 ; ii < dim ; ii++) {
    cluster_centroid[cluster_old * dim + ii] -= (X[move_point * dim + ii] - cluster_centroid[cluster_old * dim + ii]) / cluster_member_count[cluster_old];
    cluster_centroid[cluster_new * dim + ii] += (X[move_point * dim + ii] - cluster_centroid[cluster_new * dim + ii]) / cluster_member_count[cluster_new];
  }
}

void cluster_diag(int dim, int n, int k, double *X, int *cluster_assignment_index, double *cluster_centroid) {
  int cluster_member_count[MAX_CLUSTERS];

  get_cluster_member_count(n, k, cluster_assignment_index, cluster_member_count);

  printf("  Final clusters \n");
  for (int ii = 0 ; ii < k ; ii++) {
    printf(" cluster %d: members: %8d, centroid (%.1f %.1f) \n", ii, cluster_member_count[ii], cluster_centroid[ii * dim + 0], cluster_centroid[ii * dim + 1]);
  }
}

void copy_assignment_array(int n, int *src, int *tgt) {
  for (int ii = 0 ; ii < n ; ii++) {
    tgt[ii] = src[ii];
  }
}

int assignment_change_count(int n, int a[], int b[]) {
  int change_count = 0;

  for (int ii = 0 ; ii < n ; ii++) {
    if (a[ii] != b[ii]) {
      change_count++;
    }
  }

  return change_count;
}

void kmeans(
    int dim,                       // dimension of data
    double *X,                     // pointer to data
    int n,                         // number of elements
    int k,                         // number of clusters
    double *cluster_centroid,      // initial cluster centroids
    int *cluster_assignment_final  // output
) {

  double *dist = (double *) malloc(sizeof(double) * n * k);
  int *cluster_assignment_cur = (int *) malloc(sizeof(int) * n);
  int *cluster_assignment_prev = (int *) malloc(sizeof(int) * n);
  double *point_move_score = (double *) malloc(sizeof(double) * n * k);

  if (!dist || !cluster_assignment_cur || !cluster_assignment_prev || !point_move_score) {
    fail("Error allocating dist arrays");
  }

  // initial setup
  calc_all_distances(dim, n, k, X, cluster_centroid, dist);
  choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);
  copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);

  // BATCH UPDATE
  double prev_totD = BIG_double;
  int batch_iteration = 0;
  while (batch_iteration < MAX_ITERATIONS) {
//        printf("batch iteration %d \n", batch_iteration);
//        cluster_diag(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

    // update cluster centroids
    calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

    // see if we've failed to improve
    double totD = calc_total_distance(dim, n, k, X, cluster_centroid, cluster_assignment_cur);
    if (totD > prev_totD) {
      // failed to improve - currently solution worse than previous
      // restore old assignments
      copy_assignment_array(n, cluster_assignment_prev, cluster_assignment_cur);

      // recalc centroids
      calc_cluster_centroids(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

      printf("  negative progress made on this step - iteration completed (%.2f) \n", totD - prev_totD);

      // done with this phase
      break;
    }

    // save previous step
    copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_prev);

    // move all points to nearest cluster
    calc_all_distances(dim, n, k, X, cluster_centroid, dist);
    choose_all_clusters_from_distances(dim, n, k, dist, cluster_assignment_cur);

    int change_count = assignment_change_count(n, cluster_assignment_cur, cluster_assignment_prev);

    printf("%3d   %u   %9d  %16.2f %17.2f\n", batch_iteration, 1, change_count, totD, totD - prev_totD);
    fflush(stdout);

    // done with this phase if nothing has changed
    if (change_count == 0) {
      printf("  no change made on this step - iteration completed \n");
      break;
    }

    prev_totD = totD;
    batch_iteration++;
  }

  cluster_diag(dim, n, k, X, cluster_assignment_cur, cluster_centroid);

  // write to output array
  copy_assignment_array(n, cluster_assignment_cur, cluster_assignment_final);

  free(dist);
  free(cluster_assignment_cur);
  free(cluster_assignment_prev);
  free(point_move_score);
}

int main(int argc, char *argv[]) {

  FILE *inFile, *outFile;
  bool inputFile = false, outputFile = false, clusterCount = false;
  int numClusters;

  for (int j = 1 ; j < argc ; ++j) {
    if (strcmp(argv[j], "-i") == 0) {
      inputFile = true;
      if ((inFile = fopen(argv[(j + 1)], "r")) == NULL) {
        inputFile = false;
      }
    } else if (strcmp(argv[j], "-c") == 0) {
      numClusters = atoi(argv[(j + 1)]);
      clusterCount = true;
    } else if (strcmp(argv[j], "-o") == 0) {
      outputFile = true;
      outFile = fopen(argv[(j + 1)], "w");
    }
  }

  if (!inputFile) {
    printf("Incorrect input file! Using test input: \"test-kmeans.txt\" \n");
    inFile = fopen("test-kmeans.txt", "r");
  }

  if (!clusterCount) {
    printf("Cluster count not provided! Using default value = 10\n");
    numClusters = 10;
  }

  if (!outputFile) {
    printf("Output file not specified: Using default test output: \"test-kmeans-out.txt\" \n");
    outFile = fopen("test-kmeans-out.txt", "w");
  }

  int numPts = 0, i = 0;
  char ch;
  while((ch=fgetc(inFile)) != EOF) {
    if(ch == '\n')
      numPts++;
  }
  rewind(inFile);

  printf("Number of points: %d\n", numPts);
  double *points = (double *)malloc(numPts * sizeof(double));
  double *initial_centroids = (double *)malloc(numClusters * sizeof(double));
  int *assignments = (int *)malloc(numPts * sizeof(int));

  printf("Length of array: %d\n", (int)(sizeof(points)));

  char line[10];
  while (fgets(line, 10, inFile)) {
    int n = sscanf(line, "%lf", &points[i]);
    ++i;
  }

  for (int j = 0 ; j < numClusters ; ++j) {
    initial_centroids[j] = (double) j * 100.0f / numClusters;
  }

  kmeans(1, points, numPts, numClusters, initial_centroids, assignments);

  fprintf(outFile, "Data Point \t Cluster \n");
  for (int j = 0 ; j < numPts ; ++j) {
    fprintf(outFile, "%f \t %d\n", points[j], assignments[j]);
  }

  fclose(inFile);
  fclose(outFile);

  free(points);
  free(initial_centroids);
  free(assignments);

}