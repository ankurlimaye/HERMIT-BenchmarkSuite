/* Prototypes for it.c */

#include "iradoninc.h"
#include "itini.h"
#include "sparse.h"

void ReadItArgs(char *);

void InitArrays(void);
SparseMatrix *GenerateAMatrix(void);
Vector *GenerateAMatrixRow(int);
Vector *GenerateAMatrixColumn(int);

void MaskVec(Vector *, Vector *, float, float);

void SaveIteration(Vector *, int, char *);
Image *FAST_ART(SparseMatrix *, Vector *, Vector *);
Image *SLOW_ART(Vector *, Vector *);
Image *FAST_EM(SparseMatrix *, Vector *, Vector *);
Image *SLOW_EM(Vector *, Vector *);
Image *FAST_CG(SparseMatrix *, Vector *, Vector *);
Image *SLOW_CG(Vector *, Vector *);

extern itINItype itINI;





