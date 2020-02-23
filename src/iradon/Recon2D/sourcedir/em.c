/*****************************************************************
[HEADER]

This unit contains the programs associated with iterative
reconstruction with the ML-EM algorithm.

Dec. 94, JJJ and PT
*****************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "it.h"

/****************************************************************************
[NAME]
EMUpdateMultVector

[SYNOPSIS]
void EMUpdateMultVector(Vector *Xn, 
                        Vector *Nom, 
                        Vector *DeNom)

[DESCRIPTION] 
This function preforms the following EM specific update operation
$$\mathbf{a=b/c},$$
where $\mathbf a$={\tt Xn}, $\mathbf b$={\tt Nom}, $\mathbf c$={\tt DeNom}.

[USAGE]
Do not use, used internally.

[REVISION]
Jan. 95, JJJ and PT
****************************************************************************/
void EMUpdateMultVector(Vector *Xn, Vector *Nom, Vector *DeNom) {
  int n;
  float *tempX, *tempN, *tempD;

  tempX = Xn->value;
  tempN = Nom->value;
  tempD = DeNom->value;

  for (n = 0 ; n < Xn->N ; n++)
    if (fabs(tempD[n]) > 10e-6)
      tempX[n] *= tempN[n] / tempD[n];
    else
      tempX[n] = 0.0;
}

/****************************************************************************
[NAME]
FAST\_EM

[SYNOPSIS]
Image *FAST_EM(SparseMatrix *AMatrix, 
               Vector *xvector,
               Vector *bvector)


[DESCRIPTION]
This function will iterate towards a solution for the sparse system of 
equations \mb{b}=\mb{A x}, where \mb{b} is the sinogram (Radon domain)
and \mb{x} is the reconstructed image to be found. The function uses ML--EM
(Maximum Likelihood Expectation Maximization) in a fast version using a 
precalculated transformation matrix..

[USAGE]
{\tt Image=FAST\_EM(TestMatrix, TestSinogram, reference);}

Reconstructs the sinogram {\tt TestSinogram}, returns it as an image.

[REVISION]
Jan. 95, JJJ and PT
April 96 PT Support for constraints build in
****************************************************************************/
Image *FAST_EM(SparseMatrix *AMatrix, Vector *xvector, Vector *bvector) {
  int i, ARows, ACols, currentiteration, UseRefImage;
  float refxdev = 0.0, *tempXrv, *tempBv, *tempBp;
  char DiffFileName[200];
  FILE *DiffFile = NULL;
  Vector *refxvector = NULL, *bpredict, *xpredict, *ATRowSum;
  Image *Recon, *RefImage = NULL;

/*Vector *gedmask;
Image *Mask;
Mask=ReadFIF("OP_15_mask");
gedmask=ImageToVector(Mask);
FreeImage(Mask);
*/
  ARows = AMatrix->M;
  ACols = AMatrix->N;

  Print(_DNormal, "Using EM to solve %i equations with %i unknowns\n", ARows, ACols);

  UseRefImage = (strlen(itINI.RefFileName) != 0);
  if (UseRefImage != 0) {
    RefImage = ReadFIF(itINI.RefFileName);
    refxvector = ImageToVector(RefImage);
    FreeImage(RefImage);
    tempXrv = refxvector->value;
    refxdev = DeviationVector(refxvector);
    strcpy(DiffFileName, itINI.OutFileName);
    strcat(DiffFileName, ".dif");
    DiffFile = fopen(DiffFileName, "wt");
    Print(_DNormal, "Logging differences in `%s' \n", DiffFileName);
  }

  ATRowSum = SumRowSparseTMatrix(AMatrix);
  bpredict = InitVector(ARows);
  xpredict = InitVector(ACols);
  tempBv = bvector->value;

  for (currentiteration = 0 ; currentiteration < itINI.Iterations ; currentiteration++) {
    Print(_DNoLog, "Iterating: %6.2f %% done\r",
          (currentiteration + 1) * 100.0 / itINI.Iterations);
/*MaskVec(xvector,gedmask,0,10);*/
    MultSparseMatrixVector(AMatrix, xvector, bpredict);
    for (i = 0 ; i < bvector->N ; i++) {
      if (fabs(*(tempBp = &bpredict->value[i])) > 10e-6)
        *tempBp = tempBv[i] / (*tempBp);
      else
        *tempBp = 0.0;
    }

    MultSparseTMatrixVector(AMatrix, bpredict, xpredict);
    EMUpdateMultVector(xvector, xpredict, ATRowSum);

    if ((itINI.ConstrainMin >= 0) && (itINI.ConstrainMax >= 0))
      ConstrainVector(xvector, itINI.ConstrainMin, itINI.ConstrainMax);

    if (itINI.SaveIterations != 0)
      SaveIteration(xvector, currentiteration, itINI.OutFileName);

    if (UseRefImage == 1)
      fprintf(DiffFile, "%i %e\n",
              currentiteration,
              L2NormVector(refxvector, xvector, refxdev));
  }
  Print(_DNoLog, "                                                  \r");
  Recon = VectorToImage(xvector, itINI.XSamples, itINI.YSamples);

  if (UseRefImage == 1) {
    Print(_DNormal, "L2 = %9.6f \n", L2NormVector(refxvector, xvector, refxdev));
    FreeVector(refxvector);
    fclose(DiffFile);
  }
  FreeVector(bpredict);
  FreeVector(xpredict);
  FreeVector(ATRowSum);
  RenameImage(Recon, itINI.OutFileName);
  Recon->DeltaX = itINI.DeltaX;
  Recon->DeltaY = itINI.DeltaY;
  Recon->Xmin = itINI.Xmin;
  Recon->Ymin = itINI.Ymin;

  return Recon;
}

/****************************************************************************
[NAME]
SLOW\_EM

[SYNOPSIS]
Image *SLOW_EM(SparseMatrix *AMatrix, 
          Vector *xvector,
          Vector *bvector)


[DESCRIPTION]
This function will iterate towards a solution for the sparse system of 
equations \mb{b}=\mb{A x}, where \mb{b} is the sinogram (Radon domain)
and \mb{x} is the reconstructed image to be found. The function uses ML--EM
(Maximum Likelihood Expectation Maximization) in a slow version calculating
the transformaiton matrix on the fly.

[USAGE]
{\tt Image=SLOW\_EM(TestSinogram, reference);}

Reconstructs the sinogram {\tt TestSinogram}, returns it as an image.

[REVISION]
Jan. 95, JJJ and PT
April 96 Support for constraints build in PT
****************************************************************************/
Image *SLOW_EM(Vector *xvector, Vector *bvector) {
  int n, m, ARows, ACols, currentiteration, UseRefImage;
  float refxdev = 0.0, *tempXrv, *tempBv, *tempBp;
  char DiffFileName[200];
  FILE *DiffFile = NULL;
  Vector *refxvector = NULL;
  Vector *AVector, *ATVector, *RVector, *DVector;
  Vector *SVector, *YVector;
  Image *Recon, *RefImage = NULL;

  InitArrays();

  ARows = itINI.ThetaSamples * itINI.RhoSamples;
  ACols = itINI.XSamples * itINI.YSamples;
  Print(_DNormal, "Using EM to solve %i equations with %i unknowns\n",
        ARows, ACols);

  UseRefImage = (strlen(itINI.RefFileName) != 0);
  if (UseRefImage != 0) {
    RefImage = ReadFIF(itINI.RefFileName);
    refxvector = ImageToVector(RefImage);
    FreeImage(RefImage);
    tempXrv = refxvector->value;
    refxdev = DeviationVector(refxvector);
    strcpy(DiffFileName, itINI.OutFileName);
    strcat(DiffFileName, ".dif");
    DiffFile = fopen(DiffFileName, "wt");
    Print(_DNormal, "Logging differences in `%s' \n", DiffFileName);
  }

  SVector = InitVector(ACols);
  YVector = InitVector(ARows);
  DVector = InitVector(ARows);
  RVector = InitVector(ARows);
  tempBv = bvector->value;

  for (currentiteration = 0 ; currentiteration < itINI.Iterations ; currentiteration++) {
    Print(_DNoLog, "Iterating: %6.2f %% done\r",
          (currentiteration + 1) * 100.0 / itINI.Iterations);

    for (m = 0 ; m < ARows ; m++) {
      AVector = GenerateAMatrixRow(m);
      RVector->value[m] = MultVectorVector(AVector, xvector);
      FreeVector(AVector);
    }
    for (m = 0 ; m < bvector->N ; m++) {
      if (fabs(*(tempBp = &RVector->value[m])) > 10e-6)
        DVector->value[m] = tempBv[m] / (*tempBp);
      else
        DVector->value[m] = 0.0;
    }
    for (n = 0 ; n < ACols ; n++) {
      ATVector = GenerateAMatrixColumn(n);
      if (currentiteration == 0)
        for (m = 0 ; m < ARows ; m++)
          SVector->value[n] += ATVector->value[m];
      YVector->value[n] = MultVectorVector(ATVector, DVector);
      FreeVector(ATVector);
    }
    EMUpdateMultVector(xvector, YVector, SVector);

    if ((itINI.ConstrainMin >= 0) && (itINI.ConstrainMax >= 0))
      ConstrainVector(xvector, itINI.ConstrainMin, itINI.ConstrainMax);

    if (itINI.SaveIterations != 0)
      SaveIteration(xvector, currentiteration, itINI.OutFileName);

    if (UseRefImage == 1)
      fprintf(DiffFile, "%i %e\n", currentiteration,
              L2NormVector(refxvector, xvector, refxdev));
  }
  Print(_DNoLog, "                                                  \r");
  Recon = VectorToImage(xvector, itINI.XSamples, itINI.YSamples);

  if (UseRefImage == 1) {
    Print(_DNormal, "L2 = %9.6f \n", L2NormVector(refxvector, xvector, refxdev));
    FreeVector(refxvector);
    fclose(DiffFile);
  }

  FreeVector(SVector);
  FreeVector(DVector);
  FreeVector(YVector);
  FreeVector(RVector);

  RenameImage(Recon, itINI.OutFileName);
  Recon->DeltaX = itINI.DeltaX;
  Recon->DeltaY = itINI.DeltaY;
  Recon->Xmin = itINI.Xmin;
  Recon->Ymin = itINI.Ymin;

  return Recon;
}
