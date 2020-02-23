/***************************************************************
[HEADER] 

This toolbox contains functions to manipulate with, and calculate
using sparse matrices, sparse vectors and ordinary vectors.  This
library operates with the following type definitions

\begin{itemize}
\item SparseMatrix
\item SparseVector
\item Vector
\end{itemize}

all of which are defined in {\tt sparse.h}

***************************************************************/

#include "it.h"

/***************************************************************
[NAME]
InitSparseMatrix
            
[SYNOPSIS]
SparseMatrix *InitSparseMatrix(int NumRows, 
                               int NumCols)

[DESCRIPTION] 

This function initializes and allocates memory for a sparse matrix
with a specified maximum size of {\tt NumRows}$\times${\tt NumCols}.
The type {\tt SparseMatrix} is defined in {\tt sparse.h}.  All the
members of the matrix is initialized to zero.  A pointer to the new
matrix is returned.

[USAGE]
{\tt TestMatrix=InitSparseMatrix(1000,1000);}

Allocates a sparse matrix with a capacity of 1000$\times$1000 members.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
SparseMatrix *InitSparseMatrix(int NumRows, int NumCols) {
  int i;
  SparseMatrix *tempM;

  if (!(tempM = (SparseMatrix *) malloc(sizeof(SparseMatrix))))
    Error("Memory Allocation problems (InitSparseMatrix)");
  if (!(tempM->Nm = (int *) malloc(sizeof(int) * NumRows)))
    Error("Memory Allocation problems (InitSparseMatrix)");
  if (!(tempM->index = (int **) malloc(sizeof(int *) * NumRows)))
    Error("Memory Allocation problems (InitSparseMatrix)");
  if (!(tempM->value = (float **) malloc(sizeof(float *) * NumRows)))
    Error("Memory Allocation problems (InitSparseMatrix)");

  tempM->M = NumRows;
  tempM->N = NumCols;
  for (i = 0 ; i < NumRows ; i++)
    tempM->Nm[i] = 0;

  return tempM;
}

/***************************************************************
[NAME]
FreeSparseMatrix

[SYNOPSIS]
FreeSparseMatrix(SparseMatrix *MyMatrix)

[DESCRIPTION] 

This function frees the memory allocated for the sparse matrix {\tt
MyMatrix}.

[USAGE]
{\tt FreeSparseMatrix(TestMatrix);}

Frees the matrix {\tt Testmatrix}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
void FreeSparseMatrix(SparseMatrix *MyMatrix) {
  int m;

  for (m = 0 ; m < MyMatrix->M ; m++)
    if (MyMatrix->Nm[m] != 0) {
      free(MyMatrix->index[m]);
      free(MyMatrix->value[m]);
    }
  free(MyMatrix->Nm);
  free(MyMatrix->index);
  free(MyMatrix->value);
  free(MyMatrix);
}

/***************************************************************
[NAME]
InitVector

[SYNOPSIS]
Vector *InitVector(int NumElements)

[DESCRIPTION] 

This function initializes and allocates memory for a floating point
vector defined in {\tt sparse.h} as type {\tt Vector}. A pointer is
returned to the new vector.

[USAGE]
{\tt TestVector=InitVector(100);}

Initializes a normal float vector of size 100 elements. 

[REVISION]
Dec. 94, JJJ and PT\\
Oct 12 96 PT Revision if 0 elements
***************************************************************/
Vector *InitVector(int NumElements) {
  Vector *tempV;

  if (!(tempV = (Vector *) malloc(sizeof(Vector))))
    Error("Memory Allocation problems (InitVector)");
  if (NumElements > 0) {
    if (!(tempV->value = (float *) calloc(NumElements, sizeof(float))))
      Error("Memory Allocation problems (InitVector)");
  } else
    tempV->value = NULL;

  tempV->N = NumElements;

  return tempV;
}

/***************************************************************
[NAME]
FreeVector

[SYNOPSIS]
FreeVector(Vector *MyVector)

[DESCRIPTION]
Frees the memory allocated for the vector {\tt MyVector}.

[USAGE]
{\tt FreeVector(TestVector);}

Frees the vector {\tt TestVector}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
void FreeVector(Vector *MyVector) {
  if (MyVector->N > 0)
    free(MyVector->value);
  free(MyVector);
}

/***************************************************************
[NAME]
InitSparseVector

[SYNOPSIS]
SparseVector *InitSparseVector(int NumElements)

[DESCRIPTION] 

This function initializes, and allocates memory for a sparse vector
with {\tt NumElements} entries. The type {\tt SparseVector} is defined
in {\tt sparse.h}. A pointer to the new vector is returned.

[USAGE]
{\tt TestVector=SparseVector(1000);}

Allocates a sparse vector with a capacity of 1000 elements.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
SparseVector *InitSparseVector(int NumElements) {
  SparseVector *tempV;

  if (!(tempV = (SparseVector *) malloc(sizeof(SparseVector))))
    Error("Memory Allocation problems (InitSparseVector)");
  if (NumElements > 0) {
    if (!(tempV->index = (int *) malloc(sizeof(int) * NumElements)))
      Error("Memory Allocation problems (InitSparseVector)");
    if (!(tempV->value = (float *) malloc(sizeof(float) * NumElements)))
      Error("Memory Allocation problems (InitSparseVector)");
  } else {
    tempV->index = NULL;
    tempV->value = NULL;
    Print(_DDebug, "\nWarning : Zero elements initialized (InitSparseVector)\n");
  }

  tempV->N = NumElements;

  return tempV;
}

/***************************************************************
[NAME]
ConvertVector

[SYNOPSIS]
SparseVector *ConvertVector(Vector *MyVector, 
                            float MinLevel)

[DESCRIPTION] 

This function converts an ordinary vector to the sparse format,
discarding all elements with a nummeric value under {\tt MinLevel}. a
pointer to a sparse vector is returned.


[USAGE]
{\tt NewVector=ConvertVector(TestVector, 0.01);}

Converts {\tt TestVector} to the sparse vector {\tt NewVector}. All the
elements with a nummeric value under 0.01 are removed.

[REVISION] 
Dec. 94, JJJ and PT
***************************************************************/
SparseVector *ConvertVector(Vector *MyVector,
                            float MinLevel) {
  int n, tempN, tempVE, SparseElements, *ValidElements;
  float *tempF;
  SparseVector *tempV;

  tempF = MyVector->value;
  tempN = MyVector->N;
  ValidElements = IntVector(tempN);
  for (n = 0, SparseElements = 0 ; n < tempN ; n++)
    if (fabs(tempF[n]) > MinLevel) ValidElements[SparseElements++] = n;

  tempV = InitSparseVector(SparseElements);
  for (n = 0 ; n < SparseElements ; n++) {
    tempV->index[n] = (tempVE = ValidElements[n]);
    tempV->value[n] = tempF[tempVE];
  }

  free(ValidElements);
  return tempV;
}

/***************************************************************
[NAME]
InsertSparseVector

[SYNOPSIS]
void InsertSparseVector(SparseMatrix *MyMatrix, 
                        SparseVector *MyVector, 
                        int index)

[DESCRIPTION] 

This function inserts a sparse vector as a row into a sparse
matrix. This function combined with {\tt convertvector} is an easy and
convenient way of constructing sparse matrices, instead working on the
matrix itself. 

[USAGE]
{\tt InsertSparseVector(TestMatrix, TestVector, 15);}

Inserts the vector {\tt TestVector} into {\tt TestMatrix} as row 15.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
void InsertSparseVector(SparseMatrix *MyMatrix,
                        SparseVector *MyVector,
                        int index) {
  if (MyMatrix->Nm[index] != 0)
    Error("The requested row is already spoken for (InsertSparseVector)");

  MyMatrix->index[index] = MyVector->index;
  MyMatrix->value[index] = MyVector->value;
  MyMatrix->Nm[index] = MyVector->N;
  free(MyVector);
}

/***************************************************************
[NAME]
GetElement 

[SYNOPSIS]
float GetElement(SparseMatrix *MyMatrix, 
                 int m, 
                 int n)

[DESCRIPTION]

This function returns the value of the element on the {\tt m}'th row,
{\tt n}'th column in the matrix {\tt MyMatrix}.

[USAGE]
{\tt a=GetElement(TestMatrix,13,45);}

Returns the value of the elemtent in the 13'th row, 45'th column in
{\tt TestMatrix}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
float GetElement(SparseMatrix *MyMatrix, int m, int n) {

  int i;
  int *tempV, tempNm;

  tempNm = MyMatrix->Nm[m];
  tempV = MyMatrix->index[m];
  i = 0;
  while ((tempV[i] < n) && (i < tempNm)) i++;
  return (i >= MyMatrix->Nm[m]) ? (0.0) : ((n == tempV[i]) ? (MyMatrix->value[m][i]) : (0.0));
}

/***************************************************************
[NAME]
MatrixCat

[SYNOPSIS]
MatrixCat(SparseMatrix *MyM1,
          SparseMatrix *MyM2);

[DESCRIPTION] 

This function catenates the matrix {\tt MyM2} under the matrix {\tt
MyM1}, and returns the new matrix in {\tt MyM1}. The matrix {\tt MyM2}
is destroyed. The returned matrix equals $M_1=\left (
\begin{array}{c}M_1\\M_2\end{array}\right )$.

[USAGE]
{\tt MatrixCat(Test1,Test2);}

Catenates the matrix {\tt Test2} under the matrix {\tt Test1}.

[REVISION]
Dec. 94, JJJ
***************************************************************/
void MatrixCat(SparseMatrix *MyM1, SparseMatrix *MyM2) {
  SparseMatrix *NewMatrix;

  if (MyM1->N != MyM2->N)
    Error("Incompatible sizez encountered (MatrixCat)");

  NewMatrix = InitSparseMatrix((MyM1->M + MyM2->M), MyM1->N);

  memcpy(NewMatrix->value, MyM1->value, sizeof(float *) * MyM1->M);
  memcpy(NewMatrix->index, MyM1->index, sizeof(int *) * MyM1->M);
  memcpy(NewMatrix->Nm, MyM1->Nm, sizeof(int) * MyM1->M);

  memcpy(&NewMatrix->value[MyM1->M], MyM2->value, sizeof(float *) * MyM2->M);
  memcpy(&NewMatrix->index[MyM1->M], MyM2->index, sizeof(int *) * MyM2->M);
  memcpy(&NewMatrix->Nm[MyM1->M], MyM2->Nm, sizeof(int) * MyM2->M);

  free(MyM1->Nm);
  free(MyM1->index);
  free(MyM1->value);

  free(MyM2->Nm);
  free(MyM2->index);
  free(MyM2->value);
  free(MyM2);

  *MyM1 = *NewMatrix;
  free(NewMatrix);
}

/***************************************************************
[NAME]
VectorCat

[SYNOPSIS]
VectorCat(Vector *MyV1,
          Vector *MyV2);

[DESCRIPTION] 

This function catenates the vector {\tt MyV2} under the vector {\tt
MyV1}, and returns the new vector in {\tt MyV1}. The vector {\tt MyV1}
is destroyed. The returned vector equals $V_1=\left (
\begin{array}{c}V_1\\V_2\end{array}\right )$.

[USAGE]
{\tt VectorCat(Test1,Test2);}

Catenates the vector {\tt Test2} under the vector {\tt Test1}.

[REVISION]
Dec. 94, JJJ
***************************************************************/
void VectorCat(Vector *MyV1, Vector *MyV2) {
  Vector *NewVector;

  NewVector = InitVector(MyV1->N + MyV2->N);

  memcpy(NewVector->value, MyV1->value, sizeof(float) * MyV1->N);
  memcpy(&NewVector->value[MyV1->N], MyV2->value, sizeof(float) * MyV2->N);

  free(MyV1->value);
  FreeVector(MyV2);

  MyV1->value = NewVector->value;
  MyV1->N = NewVector->N;
  free(NewVector);
}

/***************************************************************
[NAME]
TransposeMatrix

[SYNOPSIS]
SparseMatrix *TransposeMatrix(SparseMatrix *MyMatrix)

[DESCRIPTION] 

This function transposes the sparse matrix {\tt MyMatrix}. A pointer
is returned to the new transposed matrix.

[USAGE]
{\tt TestMatrixT=TransposeMatrix(TestMatrix);}

Returns the transposed matrix of {\tt TestMatrix} in {\tt TestMatrixT}

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
SparseMatrix *TransposeMatrix(SparseMatrix *MyMatrix) {
  int m, n, *FillCounter, *LengthCounter, *tempI, tempNm, newm, newn;
  SparseMatrix *NewMatrix;

  float *tempV;

  FillCounter = IntVector(MyMatrix->N);
  LengthCounter = IntVector(MyMatrix->N);
  /* count the elements in every new vector! */
  for (m = 0 ; m < MyMatrix->M ; m++) {
    tempI = MyMatrix->index[m];
    tempNm = MyMatrix->Nm[m];
    for (n = 0 ; n < tempNm ; n++)
      LengthCounter[tempI[n]]++;
  }

  /* Allocate new matrix */
  NewMatrix = InitSparseMatrix(MyMatrix->N, MyMatrix->M);
  for (m = 0 ; m < NewMatrix->M ; m++)
    InsertSparseVector(NewMatrix, InitSparseVector(LengthCounter[m]), m);

  /* Insert values in new matrix */
  for (m = 0 ; m < MyMatrix->M ; m++) {
    tempI = MyMatrix->index[m];
    tempV = MyMatrix->value[m];
    tempNm = MyMatrix->Nm[m];
    for (n = 0 ; n < tempNm ; n++) {
      newn = FillCounter[(newm = tempI[n])]++;
      NewMatrix->value[newm][newn] = tempV[n];
      NewMatrix->index[newm][newn] = m;
    }
  }

  Free(FillCounter);
  Free(LengthCounter);

  return NewMatrix;
}

/***************************************************************
[NAME]
WriteSIF

[SYNOPSIS]
void WriteSIF(SparseMatrix *MyMatrix, 
              itINItype *MyitINI, 
              char* filename)

[DESCRIPTION] 

This function writes a sparse matrix to a `{\tt .sif}' file, containing
the sparse matrix values, as well as all the transformation parameters
as specified in the {\tt MyitINI} structure. This is nessesary to
ensure that the transformation matrix is only applied in the correct
images.

[USAGE]
{\tt WriteSIF(Testmatrix, TestINI, "testmatrix");}

Writes the sparse matrix {\tt TestMatrix} to the file ``testmatrix.sif''.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
void WriteSIF(SparseMatrix *MyMatrix, itINItype *MyitINI, char *filename) {
  int m;
  FILE *outfile;
  char outfilename[100];

  strcpy(outfilename, filename);
  strcat(outfilename, ".sif");

  Print(_DNormal, "WriteSIF: Writing '%s' (%dx%d) \n",
        outfilename, MyMatrix->M, MyMatrix->N);

  if (!(outfile = fopen(outfilename, "wb")))
    Error("Error opening file: '%s'", outfilename);
  fwrite(MyitINI, sizeof(itINItype), 1, outfile);
  fwrite(&MyMatrix->M, sizeof(int), 1, outfile);
  fwrite(&MyMatrix->N, sizeof(int), 1, outfile);
  fwrite(MyMatrix->Nm, sizeof(int) * MyMatrix->M, 1, outfile);

  for (m = 0 ; m < MyMatrix->M ; m++) {
    fwrite(MyMatrix->index[m], sizeof(int) * MyMatrix->Nm[m], 1, outfile);
    fwrite(MyMatrix->value[m], sizeof(float) * MyMatrix->Nm[m], 1, outfile);
  }
  fclose(outfile);
}

/***************************************************************
[NAME]
ReadSIFHeader

[SYNOPSIS]
itINItype *ReadSIFHeader(char* filename)

[DESCRIPTION] 

Reads the header from a saved matrix ({\tt .sif} file) without reading
the actual data. This function can be used to check whether the saved
matrix, complies with the actual transformation parameters, without
reading the somtimes quite large quantum of data involved (15MB is not
unusual). A pointer to a {\tt itINItype} structure is returned.

[USAGE] 

{\tt TestINI=ReadSIFHeader("LargeMatrix");}

Reads the transformation parameters from ``LargeMatrix.sif'', and
returns the parameters to {\tt TestINI}.

[REVISION]
Dec. 94, JJJ and PT\\
March 96 PT
***************************************************************/
itINItype *ReadSIFHeader(char *filename) {
  FILE *infile;
  char infilename[100];
  itINItype *MyitINI;

  MyitINI = (itINItype *) malloc(sizeof(itINItype));

  strcpy(infilename, filename);
  strcat(infilename, ".sif");
  if (!(infile = fopen(infilename, "rb")))
    Error("Error opening file: '%s'", infilename);

  fread(MyitINI, sizeof(itINItype), 1, infile);
  fclose(infile);

  return MyitINI;
}

/***************************************************************
[NAME]
ReadSIF

[SYNOPSIS]
SparseMatrix *ReadSIF(itINItype *MyitINI, 
                      char* filename)

[DESCRIPTION]

This function reads a saved sparse matrix from disk. The header is
returned in {\tt MyitINI}. A pointer to the sparse matrix is returned.

[USAGE]
{\tt TestMatrix=ReadSIF(TestINI, "testmatrix");}

Reads the file ``testmatrix.sif'' into {\tt TestMatrix}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
SparseMatrix *ReadSIF(itINItype *MyitINI, char *filename) {
  int m, N, M, *tempNm;
  FILE *infile;
  char infilename[100];
  SparseVector *tempV;
  SparseMatrix *tempM;

  strcpy(infilename, filename);
  strcat(infilename, ".sif");
  if (!(infile = fopen(infilename, "rb")))
    Error("Error opening file: '%s'", infilename);

  fread(MyitINI, sizeof(itINItype), 1, infile);
  fread(&M, sizeof(int), 1, infile);
  fread(&N, sizeof(int), 1, infile);

  tempM = InitSparseMatrix(M, N);

  Print(_DNormal, "ReadSIF: Reading '%s' (%dx%d) \n",
        infilename, tempM->M, tempM->N);

  tempNm = IntVector(tempM->M);
  fread(tempNm, sizeof(int) * tempM->M, 1, infile);

  for (m = 0 ; m < tempM->M ; m++) {
    tempV = InitSparseVector(tempNm[m]);
    fread(tempV->index, sizeof(int) * tempNm[m], 1, infile);
    fread(tempV->value, sizeof(float) * tempNm[m], 1, infile);
    InsertSparseVector(tempM, tempV, m);
  }
  fclose(infile);

  Free(tempNm);
  return (tempM);
}

/***************************************************************
[NAME]
WriteSIA

[SYNOPSIS]
void WriteSIA(SparseMatrix *MyMatrix, 
              char *filename)

[DESCRIPTION] 

This function writes a sparse matrix into an ASCII file
readable by MatLAB.

[USAGE]
{\tt WriteSIA(TestMatrix,"testmatrix");}

Writes {\tt TestMatrix} to ``testmatrix.sia''.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
void WriteSIA(SparseMatrix *MyMatrix, char *filename) {
  int m, n;
  char outfilename[100];
  FILE *outfile;

  strcpy(outfilename, filename);
  strcat(outfilename, ".sia");
  if (!(outfile = fopen(outfilename, "wt")))
    Error("Error opening file: '%s'", outfilename);

  Print(_DNormal, "WriteSIA: Writing '%s' (%dx%d) \n",
        outfilename, MyMatrix->M, MyMatrix->N);

  for (m = 0 ; m < MyMatrix->M ; m++)
    for (n = 0 ; n < MyMatrix->Nm[m] ; n++)
      fprintf(outfile, "%d %d %f\n", m + 1, MyMatrix->index[m][n] + 1, MyMatrix->value[m][n]);
  fclose(outfile);
}

/***************************************************************
[NAME]
InfoSparseMatrix

[SYNOPSIS]
void InfoSparseMatrix(SparseMatrix *MyMatrix)

[DESCRIPTION] 

This function prints out some usefull statistics about the matrix {\tt
Mymatrix}. The information includes the maximum dimensions, number of
elements and the number of empty rows.

[USAGE]
{\tt InfoSparseMatrix(TestMatrix);}

Displays information about {\tt TestMatrix}.

[REVISION]
Dec. 94, JJJ and PT\\
March 96, PT
***************************************************************/
void InfoSparseMatrix(SparseMatrix *MyMatrix) {
  unsigned int m, N, M, numelements, emptyrows;

  M = MyMatrix->M;
  N = MyMatrix->N;
  numelements = 0;
  emptyrows = 0;

  for (m = 0 ; m < MyMatrix->M ; m++) {
    if (MyMatrix->Nm[m] == 0)
      emptyrows++;
    numelements += MyMatrix->Nm[m];
  }
  Print(_DNormal, "Sparsematrix (%ux%u),\n %u of %u elements nonzero (%.2f%%),\n %u of %u rows empty.\n",
        M, N, numelements, M * N, numelements / (0.01 * N * M), emptyrows, M);
  Print(_DNormal, "Disc space required: %d bytes\n", numelements * 8 + M * 4 + 8);
}

/***************************************************************
[NAME]
SumRowSparseMatrix

[SYNOPSIS]
Vector *SumRowSparseMatrix(SparseMatrix *MyMatrix)

[DESCRIPTION] 

This function returns a vector containing the sum of each row in the
matrix {\tt MyMatrix}. The length of the vector equals the number of
rows in the matrix. A pointer to the new vector is returned.
$$sum_m=\sum_{n=0}^{N}\mathbf{A}_{m,n} \hspace{10mm} \mathrm{for}\;\; 0<m<M$$

[USAGE]
{\tt SumVector=SumRowSparseMatrix(TestMatrix);}

Returns the sum of the rows in {\tt TestMatrix} in {\tt SumVector}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
Vector *SumRowSparseMatrix(SparseMatrix *MyMatrix) {
  int m, n, tempNm;
  float *tempF, sum;
  Vector *tempV;

  tempV = InitVector(MyMatrix->M);
  for (m = 0 ; m < MyMatrix->M ; m++) {
    sum = 0.0;
    tempF = MyMatrix->value[m];
    tempNm = MyMatrix->Nm[m];
    for (n = 0 ; n < tempNm ; n++)
      sum += tempF[n];
    tempV->value[m] = sum;
  }
  return tempV;
}

/***************************************************************
[NAME]
SumRowSparseTMatrix

[SYNOPSIS]
Vector *SumRowSparseTMatrix(SparseMatrix *MyMatrix)

[DESCRIPTION] 

This function returns a vector containing the sum of each row in the
transposed matrix of {\tt MyMatrix}. (e.i. each column in {\tt
MyMatrix}. The length of the vector equals the number of columns in
the matrix. A pointer to the new vector is returned.
$$sum_m=\sum_{n=0}^{N}\mathbf{A}_{m,n}^T \hspace{10mm} \mathrm{for}\;\;
0<m<M$$

[USAGE]
{\tt SumVector=SumRowSparseTMatrix(TestMatrix);}

Returns the sum of the columns in {\tt TestMatrix} in {\tt SumVector}.

[REVISION]
Dec. 94, JJJ\\
Oct 96 PT Bug if zero vals in row\\
Oct 14 PT Bug in tempV allocation
***************************************************************/
Vector *SumRowSparseTMatrix(SparseMatrix *MyMatrix) {
  int m, n, tempNm, *tempI;
  float *tempF, *tempVv;
  Vector *tempV;

  tempV = InitVector(MyMatrix->N);
  tempVv = tempV->value;
  for (m = 0 ; m < MyMatrix->M ; m++) {
    tempF = MyMatrix->value[m];
    tempI = MyMatrix->index[m];
    tempNm = MyMatrix->Nm[m];
    for (n = 0 ; n < tempNm ; n++)
      tempVv[tempI[n]] += tempF[n];
  }
  return tempV;
}

/***************************************************************
[NAME]
SumSqRowSparseMatrix

[SYNOPSIS]
Vector *SumSqRowSparseMatrix(SparseMatrix *MyMatrix)

[DESCRIPTION]

This function returns a vector containing the square sum of each row
in the matrix {\tt MyMatrix}. The length of the vector equals the
number of rows in the matrix. A pointer to the new vector is returned.
$$sqsum_m=\sum_{n=0}^{N}\mathbf{A}_{m,n}^2 \hspace{10mm} \mathrm{for}\;\;
0<m<M$$

[USAGE]
{\tt TestVector=SumSqRowSparseMatrix(TestMatrix);}

Returns the square sum of each row in {\tt TestMatrix} in {\tt TestVector}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
Vector *SumSqRowSparseMatrix(SparseMatrix *MyMatrix) {
  int m, n, tempNm;
  float *tempF, sum;
  Vector *tempV;

  tempV = InitVector(MyMatrix->M);
  for (m = 0 ; m < MyMatrix->M ; m++) {
    sum = 0.0;
    tempF = MyMatrix->value[m];
    tempNm = MyMatrix->Nm[m];
    for (n = 0 ; n < tempNm ; n++)
      sum += tempF[n] * tempF[n];
    tempV->value[m] = sum;
  }
  return tempV;
}

/***************************************************************
[NAME]
SubtractVectors

[SYNOPSIS]
Vector *SubtractVector(Vector* vec1,Vector* vec2)

[DESCRIPTION]

This function returns a difference vector between {\tt vec1} and {\tt vec2}.

[USAGE]
{\tt TestVector=SubtractVector(vec1,vec2);}


[REVISION]
Oct 96 PT
***************************************************************/
Vector *SubtractVector(Vector *vec1, Vector *vec2) {
  int m;
  Vector *tempV;

  if ((vec1->N) != (vec1->N))
    Error("The two vectors do not match in length (SubtractVector %i %i)",
          vec1->N,
          vec2->N);

  tempV = InitVector(vec1->N);
  for (m = 0 ; m < vec1->N ; m++) {
    tempV->value[m] = (vec1->value[m]) - (vec2->value[m]);
  }
  return tempV;
}

/***************************************************************
[NAME]
MeanVector

[SYNOPSIS]
float MeanVector(Vector *MyVector);

[DESCRIPTION]

Calculates and returns the mean of the vector.

[USAGE]
{\tt f=MeanVector(Vector1);}

Calculates the mean of {\tt Vector1}.

[REVISION]
Feb. 95, JJJ
***************************************************************/
float MeanValueVector(Vector *MyVector) {
  int n;
  float *Data;
  float MV;

  Data = MyVector->value;
  for (MV = 0, n = 0 ; n < MyVector->N ; n++)
    MV += Data[n];

  MV /= (float) (MyVector->N);
  return MV;
}

/***************************************************************
[NAME]
TwoNormVector

[SYNOPSIS]
float TwoNormVector(Vector *MyVector);

[DESCRIPTION]

Calculates and returns the squared two-norm of the vector.

[USAGE]
{\tt f=TwoNormVector(Vector1);}

Calculates the two-norm of {\tt Vector1}.

[REVISION]
Oct. 96, PT
***************************************************************/
float TwoNorm(Vector *MyVector) {
  int n;
  float DV, *Data;

  Data = MyVector->value;
  for (DV = 0, n = 0 ; n < MyVector->N ; n++)
    DV += Data[n] * Data[n];

  return DV;
}

/***************************************************************
[NAME]
DeviationVector

[SYNOPSIS]
float DeviationVector(Vector *MyVector);

[DESCRIPTION]

Calculates and returns the standard deviation of the vector.

[USAGE]
{\tt f=DeviationVector(Vector1);}

Calculates the deviation of {\tt Vector1}.

[REVISION]
Feb. 95, JJJ
***************************************************************/
float DeviationVector(Vector *MyVector) {
  int n;
  float MV, DV, *Data;

  MV = MeanValueVector(MyVector);
  Data = MyVector->value;
  for (DV = 0, n = 0 ; n < MyVector->N ; n++)
    DV += (MV - Data[n]) * (MV - Data[n]);

  DV = sqrt((float) DV / MyVector->N);
  return DV;
}

/***************************************************************************
[NAME]
L2NormVector

[SYNOPSIS]
float L2NormVector(Vector *OrgVector, 
                   Vector *TestVector, 
                   float OrgDeviation)

[DESCRIPTION] 

Calculates the $L_2$--norm betveen the Original vector and the test
image. To speed up calculation, the deviation of the original image must 
be calculated first, and provided.

[USAGE]
\tc{L2=L2Norm(RefVector,TestVector,Dev);}

returns the $L_2$ norm betveen {\tt RefVector} and {\tt TestVector}.

[REVISION]
Feb. 95, JJJ
***************************************************************************/
float L2NormVector(Vector *OrgVector, Vector *TestVector, float OrgDeviation) {
  int n;
  float L2, *data1, *data2;

  if (!(OrgVector->N == TestVector->N))
    Error("Pictures must be same size (L2Norm)");

  data1 = OrgVector->value;
  data2 = TestVector->value;
  for (L2 = 0.0, n = 0 ; n < OrgVector->N ; n++)
    L2 += ((data2[n] - data1[n]) * (data2[n] - data1[n]));

  L2 = sqrt(1.0 / (OrgVector->N - 1) * L2);
  L2 /= OrgDeviation;
  return L2;
}

/***************************************************************
[NAME]
MultVectorVector

[SYNOPSIS]
float MultVectorVector(Vector *MyV1, 
                       Vector *MyV2)

[DESCRIPTION]

This funciton multiplies two ordinary vectors, and returns the
result. Appropriate checks are made to ensure the vectors are compatible.
$$f=\sum_{n=0}^Na_nb_n\hspace{10mm} \mathrm{for}\;\;0<n<N$$

[USAGE]
{\tt f=MultVectorVector(Vector1,Vector2);}

Calculates the product of the vectors {\tt Vector1} and {\tt Vector2}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
float MultVectorVector(Vector *MyV1, Vector *MyV2) {
  int n, tempN;
  float *tempV1, *tempV2, sum;

  if (MyV1->N != MyV2->N)
    Error("Incompatible sizes encountered (MultSparseVectorVector)");

  tempV1 = MyV1->value;
  tempV2 = MyV2->value;
  tempN = MyV1->N;

  sum = 0.0;
  for (n = 0 ; n < tempN ; n++)
    sum += tempV1[n] * tempV2[n];

  return sum;
}

/***************************************************************
[NAME]
MultSparseVectorVector

[SYNOPSIS]
float MultSparseVectorVector(SparseVector *MySv, 
                             Vector *MyV)

[DESCRIPTION]

This funciton multiplies a sparse vector with an ordinary vector, and
returns the result like {\tt MultVectorVector}. Appropriate checks are
made to ensure the vectors are compatible.

[USAGE]
{\tt }

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
float MultSparseVectorVector(SparseVector *MySv, Vector *MyV) {
  int n, *tempI;
  float *tempSV, *tempV, sum;

  if (MySv->N != MyV->N)
    Error("Incompatible sizes encountered (MultSparseVectorVector)");

  tempI = MySv->index;
  tempSV = MySv->value;
  tempV = MyV->value;
  sum = 0.0;
  for (n = 0 ; n < MySv->N ; n++)
    sum += tempSV[n] * tempV[tempI[n]];

  return sum;
}

/***************************************************************
[NAME]
MultSparseMatrixVector

[SYNOPSIS]
void MultSparseMatrixVector(SparseMatrix *MySm, 
                            Vector *MyV, 
                            Vector *ReturnV)

[DESCRIPTION]

This function multiplies a sparse matrix with a vector, and returns
the result in the third argument, {\tt ReturnV}. Checks are made to
ensure the three arguments are compatible in size.

[USAGE]
{\tt MultSparseMatrixVector(TestMatrix, TestVector, ReturnVector);}

Multiplies {\tt TestMatrix} with {\tt TestVector}, and returns the
result in {\tt ReturnVector}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
void MultSparseMatrixVector(SparseMatrix *MySm, Vector *MyV, Vector *ReturnV) {
  int row, n, *tempI, tempNm;
  float *tempSV, *tempV, sum;

  if ((MySm->N != MyV->N) || (MySm->M != ReturnV->N))
    Error("Incompatible sizes encountered (MultSparseMatrixVector)");

  tempV = MyV->value;
  for (row = 0 ; row < MySm->M ; row++) {
    tempI = MySm->index[row];
    tempSV = MySm->value[row];
    tempNm = MySm->Nm[row];
    sum = 0.0;
    for (n = 0 ; n < tempNm ; n++)
      sum += tempSV[n] * tempV[tempI[n]];
    ReturnV->value[row] = sum;
  }
}

/***************************************************************
[NAME]
MultSparseTMatrixVector

[SYNOPSIS]
void MultSparseTMatrixVector(SparseMatrix *MySm, 
                            Vector *MyV, 
                            Vector *ReturnV)

[DESCRIPTION]

This function multiplies the transposed of a sparse matrix with a
vector, and returns the result in the third argument, {\tt
ReturnV}. Checks are made to ensure the three arguments are compatible
in size. This function can be used to reduce storage requirements drasticly on
behalf of a slight speed reduction, because the transposed matrix
does not need to be generated and stored.

[USAGE]
{\tt MultSparseTMatrixVector(TestMatrix, TestVector, ReturnVector);}

Multiplies the transposed matrix of {\tt TestMatrix} with {\tt
TestVector}, and returns the result in {\tt ReturnVector}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
void MultSparseTMatrixVector(SparseMatrix *MySm, Vector *MyV, Vector *ReturnV) {
  int m, n, *tempI;
  float *tempV, *tempR, tempMyVv;

  if ((MySm->M != MyV->N) || (MySm->N != ReturnV->N))
    Error("Incompatible sizes encountered (MultSparseMatrixVector)");

  tempR = ReturnV->value;
  for (n = 0 ; n < ReturnV->N ; n++)
    tempR[n] = 0;

  for (m = 0 ; m < MySm->M ; m++) {
    tempI = MySm->index[m];
    tempV = MySm->value[m];
    tempMyVv = MyV->value[m];
    for (n = 0 ; n < MySm->Nm[m] ; n++)
      tempR[tempI[n]] += tempV[n] * tempMyVv;
  }
}

/***************************************************************
[NAME]
MultSparseMatrixRowVector

[SYNOPSIS]
float MultSparseMatrixRowVector(SparseMatrix *MySm, 
                                Vector *MyV, 
                                int rownr)

[DESCRIPTION]

This function multiplies a row in a sparse matrix with a vector, and
returns a float value. Essentially the same as {\tt
MultSparseVectorVector}.

[USAGE]
{\tt MultSparseMatrixRowVector(TestMatrix,TestVector,17);}

Multiplies row 17 of {\tt TestMatrix} with {\tt TestVector}.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
float MultSparseMatrixRowVector(SparseMatrix *MySm, Vector *MyV, int rownr) {
  int n, *tempI, tempNm;
  float *tempSV, *tempV, sum;

  if (MySm->N != MyV->N)
    Error("Incompatible sizes encountered (MultSparseMatrixRowVector)");

  tempI = MySm->index[rownr];
  tempSV = MySm->value[rownr];
  tempNm = MySm->Nm[rownr];
  tempV = MyV->value;
  sum = 0.0;
  for (n = 0 ; n < tempNm ; n++)
    sum += tempSV[n] * tempV[tempI[n]];
  return sum;
}

/***************************************************************
[NAME]
ConstrainVector

[SYNOPSIS]
void ConstrainVector(Vector *MyVector, 
                     float low, 
                     float high)

[DESCRIPTION]

This function constrain all the elements in a vector between two
limits. If the value is not in the range [{\tt low};{\tt high}], it is
truncated to either {\tt low} or {\tt high}.

[USAGE]
{\tt ConstrainVector(TestVector,0,1);}

Constrains alle the elements of {\tt TestVector} to the interval [0;1].

[REVISION]
Dec. 94, JJJ
***************************************************************/
void ConstrainVector(Vector *MyVector, float low, float high) {
  int n, tempN;
  float *tempV, *tempVn;

  tempN = MyVector->N;
  tempV = MyVector->value;
  for (n = 0 ; n < tempN ; n++) {
    tempVn = &tempV[n];
    if (*tempVn < low) *tempVn = low;
    else if (*tempVn > high) *tempVn = high;
  }
}

/***************************************************************
[NAME]
ImageToVector

[SYNOPSIS]
Vector *ImageToVector(Image *MyImage)

[DESCRIPTION]

This function converts an Image to a Vector. The image is of type {\tt
Image} (defined in {\tt imgtools.h}). The length of the vector is
dertermined by the images dimensions. The actual conversion is made so
that {\tt Vector->Value[m+n*M]=Image->Signal[m,n]}.

[USAGE]
{\tt TestVector=ImageToVector(TestImage);}

Converts the image {\tt TestImage} to a vector {\tt TestVector}. 

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
Vector *ImageToVector(Image *MyImage) {
  int n, m, M, N;
  float *tempV, *tempI;
  Vector *NewVector;

  if (MyImage->ArrayType == _ComplexArray)
    Error("Wrong image type (ImageToVector)");

  N = MyImage->N;
  M = MyImage->M;

  NewVector = InitVector(M * N);
  tempV = NewVector->value;
  for (m = 0 ; m < M ; m++) {
    tempI = MyImage->Signal[m];
    for (n = 0 ; n < N ; n++)
      tempV[m + n * M] = tempI[n];
  }
  return NewVector;
}

/***************************************************************
[NAME]
VectorToImage

[SYNOPSIS]
Image *VectorToImage(Vector *MyVector, 
                     int M, 
                     int N)

[DESCRIPTION]

This function converts a vector to an image of size {\tt
M}$\times${\tt N}. The sampling parameters for the image are
uninitialized. The length of the input vector must equal {\tt
M}$\times${\tt N}. The vector is converted so that {\tt
Image->Signal[m,n]=Vector->value[m+n*M]}

[USAGE]
{\tt TestImage=VectorToImage(TestVector,100,100);}

Converts the vector {\tt TestVector} to the image {\tt TestImage} of
size 100 $\times$ 100 elements.

[REVISION]
Dec. 94, JJJ and PT
***************************************************************/
Image *VectorToImage(Vector *MyVector, int M, int N) {
  int n, m;
  float *tempV, *tempI;
  Image *NewImage;

  if (MyVector->N != (M * N))
    Error("Vector size not compliant with image dimensions (VectorToImage)");

  NewImage = NewFloatImage("image", M, N, _RealArray);

  tempV = MyVector->value;
  for (m = 0 ; m < M ; m++) {
    tempI = NewImage->Signal[m];
    for (n = 0 ; n < N ; n++)
      tempI[n] = tempV[m + n * M];
  }
  return NewImage;
}

/***************************************************************
[NAME]
MeanFilterVector

[SYNOPSIS]
void MeanFilterVector(Vector *MyVector, 
                      int size,
                      int weight,
                      int width,
                      int height);

[DESCRIPTION]

This function carries out a two dimensional mean filtering of the
vector representation of an image of size {\tt width} $\times$ {\tt
height}. The mean is taken over a kernel where both side dimensions
equals {\tt size*2+1} centered over the actual pixel. The {\tt width}
and {\tt height} must be the actual size the vector is converted from,
and being converted to. The {\tt weight} factor determines a mixing of
the output vector from the original vector with the filtered vector as
$x_{new}=(1-\delta)x_{org}+\delta x_{filter}$, where $\delta$ = {\tt
weight}. The {\tt weight} must be in the interval $[0;1]$.

[USAGE]
{\tt MeanFilterVector(TestVector,2,1.0,100,100);}

Mean filters the vector {\tt TestVector}, with a kernel size of 5
$\times$ 5. The vector is a representation for a 100 $\times$ 100
image. The output is the unweighted filtered image.

[REVISION]
Dec. 94, JJJ
***************************************************************/
void MeanFilterVector(Vector *MyVector, int size, float weight, int M, int N) {
  int Kern, Kerm, n, m, nmin, nmax, mmin, mmax, count;
  float sum, area, *tempVv, *tempNv, *tempKernM;
  Vector *NewVector;

  if (MyVector->N != (M * N))
    Error("Incompatible sizez encountered (MeanFilterVector)");

  NewVector = InitVector(MyVector->N);
  area = sq(2 * size + 1);
  tempVv = MyVector->value;

  for (n = 0 ; n < N ; n++) {
    nmin = max(n - size, 0);
    nmax = min(n + size + 1, N - 1);
    tempNv = &NewVector->value[n * M];
    for (m = 0 ; m < M ; m++) {
      mmin = max(m - size, 0);
      mmax = min(m + size + 1, M - 1);
      sum = 0.0;
      for (Kern = nmin, count = 0 ; Kern < nmax ; Kern++) {
        tempKernM = &tempVv[Kern * M];
        for (Kerm = mmin ; Kerm < mmax ; Kerm++, count++)
          sum += tempKernM[Kerm];
      }
      tempNv[m] = sum / count;
    }
  }
  if (weight != 1.0) {
    tempNv = NewVector->value;
    for (n = 0 ; n < MyVector->N ; n++)
      tempNv[n] = weight * tempNv[n] + (1 - weight) * tempVv[n];
  }
  Free(MyVector->value);
  MyVector->value = NewVector->value;
  Free(NewVector);
}

/***************************************************************
[NAME]
MedianFilterVector

[SYNOPSIS]
void MedianFilterVector(Vector *MyVector, 
                        int size,
                        int median,
                        float weight,
                        int width,
			int height);

[DESCRIPTION]

This function carries out a two dimensional median filtering of the
vector representation of an image of size {\tt width} $\times$ {\tt
height}. The median is taken from a kernel where both dimensions
equals {\tt size*2+1} centered over the actual pixel. The actual value
choosen depends on {\tt median}. If this value is negative or too
large, the middle value is choosen. The {\tt width} and {\tt height}
must be the actual size the vector is converted from, and being
converted to. The {\tt weight} factor determines a mixing of the
output vector from the original vector with the filtered vector as
$x_{new}=(1-\delta)x_{org}+\delta x_{filter}$, where $\delta$ = {\tt
weight}. The {\tt weight} must be in the interval $[0;1]$.

[USAGE]
{\tt MedianFilterVector(TestVector,2,3,1.0,100,100);}

Mean filters the vector {\tt TestVector}, with a kernel size of 5
$\times$ 5. The 3'th smallest value is choosen as the filter
output. The vector is a representation of a 100 $\times$ 100
image. The output is the unweighted filtered image.

[REVISION]
Dec. 94, JJJ
***************************************************************/
void MedianFilterVector(Vector *MyVector, int size, int mediannr,
                        float weight, int M, int N) {
  int Kern, Kerm, tempKernM, n, m, nmin, nmax, mmin, mmax, count, area;
  float *data, *tempVv, *tempNv;
  Vector *NewVector;

  if (MyVector->N != (M * N))
    Error("Incompatible sizez encountered (MedianFilterVector)");

  area = sq(2 * size + 1);
  if (mediannr <= 0 || mediannr > area)
    mediannr = (area - 1) / 2;
  data = FloatVector(area);

  NewVector = InitVector(MyVector->N);
  tempVv = MyVector->value;

  for (n = 0 ; n < N ; n++) {
    nmin = max(n - size, 0);
    nmax = min(n + size + 1, N - 1);
    tempNv = &NewVector->value[n * M];
    for (m = 0 ; m < M ; m++) {
      mmin = max(m - size, 0);
      mmax = min(m + size + 1, M - 1);
      for (Kern = nmin, count = 0 ; Kern < nmax ; Kern++) {
        tempKernM = Kern * M;
        for (Kerm = mmin ; Kerm < mmax ; Kerm++)
          data[count++] = tempVv[Kerm + tempKernM];
      }
      sort(count, data);
      if (count == area)
        tempNv[m] = data[mediannr];
      else
        tempNv[m] = data[(int) (mediannr * count / area)];
    }
  }
  if (weight != 1.0) {
    tempNv = NewVector->value;
    for (n = 0 ; n < MyVector->N ; n++)
      tempNv[n] = weight * tempNv[n] + (1 - weight) * tempVv[n];
  }
  Free(data);
  Free(MyVector->value);
  MyVector->value = NewVector->value;
  Free(NewVector);
}





