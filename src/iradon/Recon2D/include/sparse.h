/* Prototypes for sparse.c */

typedef struct {
  int N;
  float *value;
} Vector;

typedef struct {
  int *index;
  float *value;
  int N;
} SparseVector;

typedef struct {
  int M;
  int N;
  int *Nm;
  int **index;
  float **value;
} SparseMatrix;

typedef struct {
  unsigned short int M;
  unsigned short int N;
  unsigned short int *Nm;
  unsigned short int **index;
  unsigned char **value;
} SmallSparseMatrix;

typedef struct {
  unsigned short int *index;
  unsigned char *value;
  unsigned short int N;
} SmallSparseVector;

SparseMatrix *InitSparseMatrix(int, int);

SparseVector *InitSparseVector(int);
void FreeSparseMatrix(SparseMatrix *);

Vector *InitVector(int);
void FreeVector(Vector *);

void InsertSparseVector(SparseMatrix *, SparseVector *, int);
SparseVector *ConvertVector(Vector *, float);

SparseMatrix *TransposeMatrix(SparseMatrix *);
float GetElement(SparseMatrix *, int, int);

float NormSparseVector(SparseVector *);
float NormSparseMatrix(SparseMatrix);

float MeanValueVector(Vector *);
float DeviationVector(Vector *);
float TwoNorm(Vector *);

float L2NormVector(Vector *, Vector *, float);

void PrintSparseMatrix(SparseMatrix *);
void InfoSparseMatrix(SparseMatrix *);

Image *VectorToImage(Vector *, int, int);
Vector *ImageToVector(Image *);
void WriteSIA(SparseMatrix *, char *);
void WriteSIF(SparseMatrix *, itINItype *, char *);
SparseMatrix *ReadSIF(itINItype *, char *);
itINItype *ReadSIFHeader(char *);

float MultVectorVector(Vector *, Vector *);
float MultSparseVectorVector(SparseVector *, Vector *);
float MultSparseMatrixRowVector(SparseMatrix *, Vector *, int);
void MultSparseMatrixVector(SparseMatrix *, Vector *, Vector *);
void MultSparseTMatrixVector(SparseMatrix *, Vector *, Vector *);

Vector *SumRowSparseMatrix(SparseMatrix *);
Vector *SumRowSparseTMatrix(SparseMatrix *);
Vector *SumSqRowSparseMatrix(SparseMatrix *);

void ConstrainVector(Vector *, float, float);
void MeanFilterVector(Vector *, int, float, int, int);
void MedianFilterVector(Vector *, int, int, float, int, int);

void MatrixCat(SparseMatrix *, SparseMatrix *);
void VectorCat(Vector *, Vector *);

Vector *SubtractVector(Vector *, Vector *);

void RegulateL1Matrix(SparseMatrix *);
void RegulateMatrix(SparseMatrix *);







