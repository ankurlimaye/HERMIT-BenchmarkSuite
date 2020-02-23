/*****************************************************************
[HEADER]

This unit contains the programs associated with generation of the
transformation matrix, in iterative reconstruction algorithms.

Dec. 94, JJJ and PT\\
March 96, PT
*****************************************************************/

#include "it.h"
static int XSamples, YSamples, RhoSamples, ThetaSamples;
static float Xmin, Ymin, DeltaX, DeltaY, RhoMin, ThetaMin, DeltaRho, DeltaTheta;
static float *xar, *yar, *rhoar, *costhetaar, *sinthetaar, *mintrigar, *thetaar;
extern char *RadonKernelar;

/*****************************************************************
[NAME]
initArrays

[SYNOPSIS]
void InitArrays(void)

[DESCRIPTION]
This function will initialize several arrays to help speed up the
Radon matrix generation.

[USAGE]

Do not use, used internally.

[REVISION]
Dec. 94, JJJ and PT
*****************************************************************/
void InitArrays(void) {
  int m, n, r, t;
  float theta, ct, st, absct, absst;

  XSamples = itINI.XSamples;
  YSamples = itINI.YSamples;
  ThetaSamples = itINI.ThetaSamples;
  RhoSamples = itINI.RhoSamples;
  Xmin = itINI.Xmin;
  Ymin = itINI.Ymin;
  RhoMin = itINI.RhoMin;
  ThetaMin = itINI.ThetaMin;
  DeltaX = itINI.DeltaX;
  DeltaY = itINI.DeltaY;
  DeltaRho = itINI.DeltaRho;
  DeltaTheta = itINI.DeltaTheta;

  xar = FloatVector(XSamples);
  yar = FloatVector(YSamples);
  rhoar = FloatVector(RhoSamples);
  thetaar = FloatVector(ThetaSamples);
  costhetaar = FloatVector(ThetaSamples);
  sinthetaar = FloatVector(ThetaSamples);
  mintrigar = FloatVector(ThetaSamples);

  for (t = 0 ; t < itINI.ThetaSamples ; t++) {
    theta = t * itINI.DeltaTheta;
    thetaar[t] = theta;
    costhetaar[t] = ct = cos(theta);
    sinthetaar[t] = st = sin(theta);
    absct = fabs(ct);
    absst = fabs(st);
    mintrigar[t] = (absct > absst) ? (1 / absct) : (1 / absst);
  }
  for (m = 0 ; m < itINI.XSamples ; m++)
    xar[m] = itINI.Xmin + m * itINI.DeltaX;
  for (n = 0 ; n < itINI.YSamples ; n++)
    yar[n] = itINI.Ymin + n * itINI.DeltaX;
  for (r = 0 ; r < itINI.RhoSamples ; r++)
    rhoar[r] = itINI.RhoMin + r * itINI.DeltaRho;
}

/*******************************************************************
[NAME]
FreeArrays

[SYNOPSIS]
void FreeArrays(void)

[DESCRIPTION]
Frees the internal help arrays, called by {\tc GenerateAMAtrix}.

[USAGE]
Do not use, used internally.

[REVISION]
Dec. 94 PT, and JJJ
*******************************************************************/
void FreeArrays(void) {
  free(xar);
  free(yar);
  free(rhoar);
  free(thetaar);
  free(costhetaar);
  free(sinthetaar);
  free(mintrigar);
}

/***********************************************************************
[NAME]
ReadAMatrix

[SYNOPSIS]
SparseMatrix *ReadAMatrix()

[DESCRIPTION] 

This function is called before the actual calculation of th
A-Matrix. If an old matrix is found on the disc, it is checked whether
the transformation parameters used to generate it equals the actual
parameters. If so, the matrix is read, and used instead of calculating
it again.

[USAGE]
Do not use, used internally. 

[REVISION] 
Dec. 94, JJJ\\ 
March 95, PT Now writes which parameter made
             the system matrix be computed.
***********************************************************************/
SparseMatrix *ReadAMatrix(void) {
  char tempFN[200];
  itINItype *tempINI;
  SparseMatrix *NewMatrix;

  strcpy(tempFN, itINI.KernelFileName);
  NewMatrix = NULL;

  if ((fopen(strcat(tempFN, ".sif"), "rb"))) {
    Print(_DNormal, "Transformation matrix found...\n");
    tempINI = ReadSIFHeader(itINI.KernelFileName);
    if (memcmp(tempINI, &itINI, sizeof(int) * 15) == 0) {
      Print(_DNormal, "Using `%s' as transformation matrix [method %i]\n",
            itINI.KernelFileName, itINI.RadonKernel);
      NewMatrix = ReadSIF(tempINI, itINI.KernelFileName);
    } else {
      NewMatrix = NULL;

      if (tempINI->XSamples != itINI.XSamples)
        Print(_DNormal, "XSamples differs (%i <-> %i)\n",
              tempINI->XSamples, itINI.XSamples);
      if (tempINI->YSamples != itINI.YSamples)
        Print(_DNormal, "YSamples differs (%i <-> %i)\n",
              tempINI->YSamples, itINI.YSamples);
      if (tempINI->DeltaX != itINI.DeltaX)
        Print(_DNormal, "DeltaX differs (%f <-> %f)\n",
              tempINI->DeltaX, itINI.DeltaX);
      if (tempINI->DeltaY != itINI.DeltaY)
        Print(_DNormal, "DeltaY differs (%f <-> %f)\n",
              tempINI->DeltaY, itINI.DeltaY);
      if (tempINI->Xmin != itINI.Xmin)
        Print(_DNormal, "Xmin differs (%f <-> %f)\n",
              tempINI->Xmin, itINI.Xmin);
      if (tempINI->Ymin != itINI.Ymin)
        Print(_DNormal, "Ymin differs (%f <-> %f)\n",
              tempINI->Ymin, itINI.Ymin);
      if (tempINI->ThetaSamples != itINI.ThetaSamples)
        Print(_DNormal, "ThetaSamples differs (%i <-> %i)\n",
              tempINI->ThetaSamples, itINI.ThetaSamples);
      if (tempINI->RhoSamples != itINI.RhoSamples)
        Print(_DNormal, "RhoSamples differs (%i <-> %i)\n",
              tempINI->RhoSamples, itINI.RhoSamples);
      if (tempINI->DeltaRho != itINI.DeltaRho)
        Print(_DNormal, "DeltaRho differs (%f <-> %f)\n",
              tempINI->DeltaRho, itINI.DeltaRho);
      if (tempINI->DeltaTheta != itINI.DeltaTheta)
        Print(_DNormal, "DeltaTheta differs (%f <-> %f)\n",
              tempINI->DeltaTheta, itINI.DeltaTheta);
      if (tempINI->ThetaMin != itINI.ThetaMin)
        Print(_DNormal, "ThetaMin differs (%f <-> %f)\n",
              tempINI->ThetaMin, itINI.ThetaMin);
      if (tempINI->RhoMin != itINI.RhoMin)
        Print(_DNormal, "RhoMin differs (%f <-> %f)\n",
              tempINI->RhoMin, itINI.RhoMin);
      if (tempINI->RadonKernel != itINI.RadonKernel)
        Print(_DNormal, "RadonKernel differs (%i <-> %i)\n",
              tempINI->RadonKernel, itINI.RadonKernel);
      if (tempINI->LowestALevel != itINI.LowestALevel)
        Print(_DNormal, "LowestALevel differs (%f <-> %f)\n",
              tempINI->LowestALevel, itINI.LowestALevel);
      if (tempINI->OverSamp != itINI.OverSamp)
        Print(_DNormal, "OverSamp differs (%i <-> %i)\n",
              tempINI->OverSamp, itINI.OverSamp);
    }
    Free(tempINI);
  }
  if (!(NewMatrix))
    Print(_DNormal,
          "Failed to use old matrix, calculating new matrix [method %i]\n",
          itINI.RadonKernel);
  return NewMatrix;
}

/***********************************************************************
[NAME]
RegulateMatrix

[SYNOPSIS]
void RegulateMatrix(SparseMatrix *AMatrix);

[DESCRIPTION] 

This function appends a regularisation matrix to the matrix {\tt
AMatrix}. The regularisation is based on the ordinary laplace operator
-1,-1,4,-1,-4 (as array) times the {\tt Regularisation} parameter.

[USAGE]
{\tt RegulateMatrix(TestMatrix);}

Calculates and appends a regularisation matrix to {\tt TestMatrix}.

[REVISION]
Jan. 95, JJJ\\
April 7 PT Includes regularization parameter now
***********************************************************************/
void RegulateMatrix(SparseMatrix *AMatrix) {
  int x, y, i;
  float Laplace[5] = {-1.0, -1.0, 4.0, -1.0, -1.0};
  int *tempVi;
  float *tempVv;
  SparseMatrix *RMatrix;
  SparseVector *vec;

  Print(_DDebug, "RegulateMatrix\n");

  XSamples = itINI.XSamples;
  YSamples = itINI.YSamples;
  for (i = 0 ; i < 5 ; i++)
    Laplace[i] *= itINI.Regularization;
  RMatrix = InitSparseMatrix(XSamples * YSamples, XSamples * YSamples);
  for (y = 1 ; y < (YSamples - 1) ; y++) {
    for (x = 1 ; x < (XSamples - 1) ; x++) {
      vec = InitSparseVector(5);
      tempVi = vec->index;
      tempVv = vec->value;
      tempVi[0] = (y - 1) * XSamples + x;
      tempVi[1] = y * XSamples + x - 1;
      tempVi[2] = y * XSamples + x;
      tempVi[3] = y * XSamples + x + 1;
      tempVi[4] = (y + 1) * XSamples + x;
      memcpy(vec->value, &Laplace, sizeof(float) * 5);
      InsertSparseVector(RMatrix, vec, x + y * XSamples);
    }
  }
  Print(_DNormal, "Adding regulating matrix to A-matrix\n");
  MatrixCat(AMatrix, RMatrix);
}

/***********************************************************************
[NAME]
RegulateL1Matrix

[SYNOPSIS]
void RegulateL1Matrix(SparseMatrix *AMatrix);

[DESCRIPTION] 

This function appends a regularisation matrix to the matrix {\tt
AMatrix}. The regularisation is based on the identity matrix
the {\tt Regularisation} parameter.

[USAGE]
{\tt RegulateL1Matrix(TestMatrix);}

Calculates and appends a regularisation matrix to {\tt TestMatrix}.

[REVISION]
Oct. 96 PT
***********************************************************************/
void RegulateL1Matrix(SparseMatrix *AMatrix) {
  int x, y;
  int *tempVi;
  float *tempVv;
  float Reg[1] = {1.0};
  SparseMatrix *RMatrix;
  SparseVector *vec;

  Print(_DDebug, "RegulateMatrix\n");
  Print(_DNormal, "Appending Regularization term\n");

  XSamples = itINI.XSamples;
  YSamples = itINI.YSamples;
  Reg[0] *= itINI.Regularization;
  RMatrix = InitSparseMatrix(XSamples * YSamples, XSamples * YSamples);
  for (y = 1 ; y < YSamples ; y++) {
    for (x = 1 ; x < XSamples ; x++) {
      vec = InitSparseVector(1);
      tempVi = vec->index;
      tempVv = vec->value;
      tempVi[0] = y * XSamples + x;
      memcpy(vec->value, &Reg, sizeof(float *));
      InsertSparseVector(RMatrix, vec, x + y * XSamples);
    }
  }
  Print(_DNormal, "Adding regulating matrix to A-matrix\n");
  MatrixCat(AMatrix, RMatrix);
}

/***********************************************************************
[NAME]
SquareRadon

[SYNOPSIS]
float SquareRadon(float rho,
                  float theta,
                  float x, 
                  float y)

[DESCRIPTION]

This function calculates the radon transform for a given line (defined
by $\rho$ and $\theta$), through the pixel in ({\tt x,y}). Scaling is
done accordingly to the actual transformation parameters. This
function is closely related to radon transformation of a square, as
found in {\tt radonana.c}. 
Now the theta is transferred to the function.

[USAGE]
Do not use, used internally.

[REVISION]
Dec. 94, JJJ and PT
Oct. 11 Now transfer theta - Not index.
***********************************************************************/
float SquareRadon(float rho, float theta, float x, float y) {
  float base, tantheta, xp1, yp1, xm1, sintheta, costheta;

  costheta = cos(theta);
  sintheta = sin(theta);
  rho = (rho - (x * costheta + y * sintheta)) / DeltaX * 2;

  if (fabs(rho) > sqrt((float) 2))
    base = 0.0;
  else {
    if (rho < 0) {
      rho = -rho;
      theta += PI;
    }
    if (theta >= PI) theta = 2 * PI - theta;
    if (theta >= PI * 0.5) theta = PI - theta;
    if (theta >= PI * 0.25) theta = PI * 0.5 - theta; /* Now theta <PI/4 */
    costheta = cos(theta);
    sintheta = sin(theta);
    tantheta = sintheta / costheta;
    xp1 = rho / costheta - tantheta;
    xm1 = rho / costheta + tantheta;
    if (theta < 1e-20)
      if (rho < 1)
        base = DeltaX;
      else
        base = 0.0;
    else {
      yp1 = (rho - costheta) / sintheta;
      if (xp1 >= 1)
        base = 0.0;
      else if (xm1 <= 1)
        base = DeltaX / costheta;
      else
        base = hypot((1 - xp1), (1 - yp1)) * DeltaX / 2;
    }
  }
  return base;
}

/***********************************************************************
[NAME]
GenerateAMatrix

[SYNOPSIS]
SparseMatrix *GenerateAMatrix()

[DESCRIPTION] 

This function generates the transformation matrix (AMtrix) used by the
iterative reconstruction methods. Uses the parameters in the {\tt
itINI} structure for the generation.

[USAGE]
{\tt Matrix = GenerateAMatrix();}

Initializes the matrix {\tt TestMatrix} as a transformation matrix.

[REVISION]
Dec. 94, JJJ and PT\\
March 96 PT (Corrected offset in kernels)\\
April 2, 96 PT Uses itini.h defs
***********************************************************************/
SparseMatrix *GenerateAMatrix() {
  int orho, otheta, m, n, r, t, intn, intm;
  float psi, pidivDeltaX, rho, mi, costheta, sintheta, of;
  float cst1, cst2, cst3, realm, realn, deltan, deltam, onedivsq2, theta;
  float *tempV;
  float normfactor, OverSampDeltaTheta, OverSampDeltaRho, absct, absst;
  Vector *vec;
  SparseMatrix *NewMatrix;
  normfactor = 1.0 / ((1 + 2 * (float) itINI.OverSamp) * (1 + 2 * itINI.OverSamp));
  pidivDeltaX = PI / DeltaX;
  onedivsq2 = 1.0 / sqrt(2.0);

  if (!(NewMatrix = ReadAMatrix())) {
    InitArrays();
    NewMatrix = InitSparseMatrix(RhoSamples * ThetaSamples, XSamples * YSamples);
    OverSampDeltaTheta = itINI.DeltaTheta / (1 + 2 * itINI.OverSamp);
    OverSampDeltaRho = itINI.DeltaRho / (1 + 2 * itINI.OverSamp);

    for (t = 0 ; t < ThetaSamples ; t++) {
      Print(_DNoLog, "Calculating: %4d of %d \r", t, ThetaSamples);
      for (r = 0 ; r < RhoSamples ; r++) {
        vec = InitVector(NewMatrix->N);
        tempV = vec->value;
        for (otheta = -itINI.OverSamp ; otheta <= itINI.OverSamp ; otheta++) {
          theta = thetaar[t] + otheta * OverSampDeltaTheta;
          costheta = cos(theta);
          sintheta = sin(theta);
          absct = fabs(costheta);
          absst = fabs(sintheta);
          mi = (absct > absst) ? (1 / absct) : (1 / absst);

          if (absst > onedivsq2) {
            cst1 = Ymin * sintheta + Xmin * costheta;
            cst2 = 1.0 / (sintheta * DeltaY);
            cst3 = DeltaX / DeltaY * costheta / sintheta;
          } else {
            cst1 = Xmin * costheta + Ymin * sintheta;
            cst2 = 1.0 / (costheta * DeltaX);
            cst3 = DeltaY / DeltaX * sintheta / costheta;
          }

          for (orho = -itINI.OverSamp ; orho <= itINI.OverSamp ; orho++) {
            rho = rhoar[r] + orho * OverSampDeltaRho;
            switch (itINI.RadonKernel) {
              case _SINC : /* Low Pass Interpolation */
                for (m = 0 ; m < XSamples ; m++) {
                  of = rho - xar[m] * costheta;
                  for (n = 0 ; n < YSamples ; n++) {
                    psi = pidivDeltaX * (of - yar[n] * sintheta);
                    if (psi < 0) psi = -psi;
                    if (psi > 1e-5)
                      tempV[m + n * XSamples] += DeltaX * sin(psi * mi) / psi;
                    else
                      tempV[m + n * XSamples] += DeltaX;
                  }
                }
                break;
              case _NN : /* Nearest neighbour */
                for (m = 0 ; m < XSamples ; m++) {
                  of = rho - xar[m] * costheta;
                  for (n = 0 ; n < YSamples ; n++) {
                    psi = of - yar[n] * sintheta;
                    if ((2 * fabs(psi * costheta) < DeltaX) && (2 * fabs(psi * sintheta) < DeltaX))
                      tempV[m + n * XSamples] += DeltaX;
                  }
                }
                break;
              case _RNN: /* ray driven nearest neighbor */
                if (absst > onedivsq2) {
                  realn = (rho - cst1) * cst2;
                  for (m = 0 ; m < XSamples ; m++) {
                    intn = (int) floor(realn + 0.5);
                    if ((intn >= 0) && (intn < YSamples))
                      tempV[m + intn * XSamples] += DeltaX;
                    realn -= cst3;
                  }
                } else {
                  realm = (rho - cst1) * cst2;
                  for (n = 0 ; n < YSamples ; n++) {
                    intm = (int) floor(realm + 0.5);
                    if ((intm >= 0) && (intm < XSamples))
                      tempV[intm + n * XSamples] += DeltaX;
                    realm -= cst3;
                  }
                }
                break;
              case _RL: /* ray driven linear interpolation */
                if (absst > onedivsq2) {
                  realn = (rho - cst1) * cst2;
                  for (m = 0 ; m < XSamples ; m++) {
                    intn = (int) floor(realn);
                    if ((intn >= 0) && (intn < YSamples - 1)) {
                      deltan = realn - intn;
                      tempV[m + intn * XSamples] += DeltaX * (1 - deltan);
                      tempV[m + (intn + 1) * XSamples] += DeltaX * deltan;
                    }
                    realn -= cst3;
                  }
                } else {
                  realm = (rho - cst1) * cst2;
                  for (n = 0 ; n < YSamples ; n++) {
                    intm = (int) floor(realm);
                    if ((intm >= 0) && (intm < XSamples - 1)) {
                      deltam = realm - intm;
                      tempV[intm + n * XSamples] += DeltaX * (1 - deltam);
                      tempV[(intm + 1) + n * XSamples] += DeltaX * deltam;
                    }
                    realm -= cst3;
                  }
                }
                break;
              case _P1:        /* pixel driven interpolation via. lines */
                if (sintheta > onedivsq2) {
                  realn = (-Ymin + (rho - Xmin * costheta) / sintheta) / DeltaY;
                  for (m = 0 ; m < XSamples ; m++) {
                    intn = (int) floor(realn);
                    if ((intn >= 0) && (intn < (YSamples - 1))) {
                      tempV[m + intn * XSamples] += SquareRadon(rho, theta, xar[m], yar[intn]);
                      tempV[m + (intn + 1) * XSamples] += SquareRadon(rho, theta, xar[m], yar[intn + 1]);
                    }
                    realn -= cst3;
                  }
                } else {
                  realm = (-Xmin + (rho - Ymin * sintheta) / costheta) / DeltaX;
                  for (n = 0 ; n < YSamples ; n++) {
                    intm = (int) floor(realm);
                    if ((intm >= 0) && (intm < (XSamples - 1))) {
                      tempV[intm + n * XSamples] += SquareRadon(rho, theta, xar[intm], yar[n]);
                      tempV[intm + 1 + n * XSamples] += SquareRadon(rho, theta, xar[intm + 1], yar[n]);
                    }
                    realm -= cst3;
                  }
                }
                break;
              case _P2:        /* pixel driven interpolation for all points */
                for (m = 0 ; m < XSamples ; m++)
                  for (n = 0 ; n < YSamples ; n++)
                    tempV[m + n * XSamples] += SquareRadon(rho, theta, xar[m], yar[n]);
                break;
              default :Error("Wrong interpolation type detected");
            }
          } /* orho */
        } /* otheta */
        for (m = 0 ; m < XSamples * YSamples ; m++)
          tempV[m] *= normfactor;

        InsertSparseVector(NewMatrix,
                           ConvertVector(vec, DeltaX * itINI.LowestALevel),
                           t + r * ThetaSamples);
/*        printf("\n%i %i %i %i\n",
               t+r*ThetaSamples,
               NewMatrix->Nm[t+r*ThetaSamples],
               NewMatrix->index[t+r*ThetaSamples],
               NewMatrix->value[t+r*ThetaSamples]);*/
        FreeVector(vec);
      } /* rho */
    } /* theta */

    FreeArrays();
    Print(_DNoLog, "                                               \r");
    InfoSparseMatrix(NewMatrix);
    WriteSIF(NewMatrix, &itINI, itINI.KernelFileName);
  }

  if (itINI.SaveMatLab) WriteSIA(NewMatrix, itINI.KernelFileName);
  if (itINI.Regularization > 1e-9)
    RegulateMatrix(NewMatrix);
  else if (itINI.Regularization < -1e-9) {
    itINI.Regularization *= -1;
    RegulateL1Matrix(NewMatrix);
  }
  InfoSparseMatrix(NewMatrix);
  return NewMatrix;
}

/***********************************************************************
[NAME]
GenerateAMatrixRow

[SYNOPSIS]
Vector *GenerateAMatrixRow(int rowindex)

[DESCRIPTION] 

This function generates a row of the the transformation matrix (AMtrix) 
used by the iterative reconstruction methods. Uses the parameters 
in the {\tt itINI} structure for the generation.

[USAGE]
{\tt Vector = GenerateAMatrixRow(1);}

Initializes the vector {\tt Vector} as row 1 in a transformation matrix.

[REVISION]
March 96, JJJ and PT\\
April 2, 96 PT Uses itini.h defs
***********************************************************************/
Vector *GenerateAMatrixRow(int rowindex) {
  int m, n, r, t, intn, intm, N;
  float cst3, psi, pidivDeltaX, rho, mi;
  float onedivsq2, costheta, sintheta, of, realm, realn, deltan, deltam;
  Vector *NewVector;
  float *tempV;

  N = XSamples * YSamples;

  /* init vector to zero */
  NewVector = InitVector(N);

  tempV = NewVector->value;
  pidivDeltaX = PI / DeltaX;
  onedivsq2 = 1.0 / sqrt(2.0);

  t = rowindex % ThetaSamples;
  r = rowindex / ThetaSamples;

  sintheta = sinthetaar[t];
  costheta = costhetaar[t];
  mi = mintrigar[t];
  rho = rhoar[r];

  if (sintheta > onedivsq2)
    cst3 = costheta / sintheta;
  else
    cst3 = sintheta / costheta;

  switch (itINI.RadonKernel) {
    case _SINC: /* Low Pass Interpolation */
      for (m = 0 ; m < XSamples ; m++) {
        of = rho - xar[m] * costheta;
        for (n = 0 ; n < YSamples ; n++) {
          psi = pidivDeltaX * (of - yar[n] * sintheta);
          if (psi < 0) psi = -psi;
          if (psi > 1e-5)
            tempV[m + n * XSamples] = DeltaX * sin(psi * mi) / psi;
          else
            tempV[m + n * XSamples] = DeltaX;
        }
      }
      break;
    case _NN: /* Nearest neighbor */
      for (m = 0 ; m < XSamples ; m++) {
        of = rho - xar[m] * costheta;
        for (n = 0 ; n < YSamples ; n++) {
          psi = of - yar[n] * sintheta;
          if ((2 * fabs(psi * costheta) < DeltaX) && (2 * fabs(psi * sintheta) < DeltaX))
            tempV[m + n * XSamples] = DeltaX;
          else
            tempV[m + n * XSamples] = 0.0;
        }
      }
      break;
    case _RNN: /* ray driven nearest neighbor */
      if (sintheta > onedivsq2) {
        realn = (-Ymin + (rho - Xmin * costheta) / sintheta) / DeltaY;
        for (m = 0 ; m < XSamples ; m++) {
          intn = (int) floor(realn + 0.5);
          if ((intn >= 0) && (intn < YSamples))
            tempV[m + intn * XSamples] = DeltaX;
          realn -= cst3;
        }
      } else {
        realm = (-Xmin + (rho - Ymin * sintheta) / costheta) / DeltaX;
        for (n = 0 ; n < YSamples ; n++) {
          intm = (int) floor(realm + 0.5);
          if ((intm >= 0) && (intm < XSamples))
            tempV[intm + n * XSamples] = DeltaX;
          realm -= cst3;
        }
      }
      break;
    case _RL: /* ray driven linear interpolation */
      if (sintheta > onedivsq2) {
        realn = (-Ymin + (rho - Xmin * costheta) / sintheta) / DeltaY;
        for (m = 0 ; m < XSamples ; m++) {
          intn = (int) floor(realn);
          if ((intn >= 0) && (intn < YSamples - 1)) {
            tempV[m + intn * XSamples] = DeltaX * (1 - (deltan = realn - intn));
            tempV[m + (intn + 1) * XSamples] = DeltaX * deltan;
          }
          realn -= cst3;
        }
      } else {
        realm = (-Xmin + (rho - Ymin * sintheta) / costheta) / DeltaX;
        for (n = 0 ; n < YSamples ; n++) {
          intm = (int) floor(realm);
          if ((intm >= 0) && (intm < XSamples - 1)) {
            tempV[intm + n * XSamples] = DeltaX * (1 - (deltam = realm - intm));
            tempV[(intm + 1) + n * XSamples] = DeltaX * deltam;
          }
          realm -= cst3;
        }
      }
      break;
    case _P1:        /* pixel driven interpolation via. lines */
      if (sintheta > onedivsq2) {
        realn = (-Ymin + (rho - Xmin * costheta) / sintheta) / DeltaY;
        for (m = 0 ; m < XSamples ; m++) {
          intn = (int) floor(realn);
          if ((intn >= 0) && (intn < (YSamples - 1))) {
            tempV[m + intn * XSamples] = SquareRadon(rho, thetaar[t], xar[m], yar[intn]);
            tempV[m + (intn + 1) * XSamples] = SquareRadon(rho, thetaar[t], xar[m], yar[intn + 1]);
          }
          realn -= cst3;
        }
      } else {
        realm = (-Xmin + (rho - Ymin * sintheta) / costheta) / DeltaX;
        for (n = 0 ; n < YSamples ; n++) {
          intm = (int) floor(realm);
          if ((intm >= 0) && (intm < (XSamples - 1))) {
            tempV[intm + n * XSamples] = SquareRadon(rho, thetaar[t], xar[intm], yar[n]);
            tempV[intm + 1 + n * XSamples] = SquareRadon(rho, thetaar[r], xar[intm + 1], yar[n]);
          }
          realm -= cst3;
        }
      }
      break;
    case _P2:        /* pixel driven interpolation for all points */
      for (m = 0 ; m < XSamples ; m++)
        for (n = 0 ; n < YSamples ; n++)
          tempV[m + n * XSamples] = SquareRadon(rho, thetaar[t], xar[m], yar[n]);
      break;
    default :Error("Wrong interpolation type detected");
  }
  return NewVector;
}

/***********************************************************************
[NAME]
GenerateAMatrixColumn

[SYNOPSIS]
Vector *GenerateAMatrixColumn(int colindex)

[DESCRIPTION] 

This function generates a column of the the transformation matrix (AMatrix) 
used by the iterative reconstruction methods. Uses the parameters 
in the {\tt itINI} structure for the generation.

[USAGE]
{\tt Vector = GenerateAMatrixRow(123);}

Initializes the vector {\tt Vector} as column 123 in a transformation matrix.

[REVISION]
March 96, JJJ and PT\\
April 2, 96 PT Corrected bug in Ray driven linear\\
April 2, 96 PT Uses itini.h defs
***********************************************************************/
Vector *GenerateAMatrixColumn(int colindex) {
  int m, n, r, t, intn, intm, M;
  float cst1, cst2, psi, pidivDeltaX, costheta, sintheta, of;
  float onedivsq2, realm, realn, deltan, deltam;
  float *tempV;
  Vector *NewVector;

  m = colindex % XSamples;
  n = colindex / XSamples;
  M = RhoSamples * ThetaSamples;
  NewVector = InitVector(M);
  tempV = NewVector->value;

  pidivDeltaX = PI / DeltaX;
  onedivsq2 = 1.0 / sqrt(2.0);

  switch (itINI.RadonKernel) {
    case _SINC: /* Low Pass Interpolation */
      for (t = 0 ; t < ThetaSamples ; t++) {
        for (r = 0 ; r < RhoSamples ; r++) {
          of = rhoar[r] - xar[m] * costhetaar[t];
          psi = pidivDeltaX * (of - yar[n] * sinthetaar[0]);
          if (psi < 0) psi = -psi;
          if (psi > 1e-5)
            tempV[t + r * ThetaSamples] = DeltaX * sin(psi * mintrigar[t]) / psi;
          else
            tempV[t + r * ThetaSamples] = DeltaX;
        }
      }
      break;

    case _NN: /* Nearest neighbor */
      for (t = 0 ; t < ThetaSamples ; t++) {
        sintheta = sinthetaar[t];
        costheta = costhetaar[t];
        for (r = 0 ; r < RhoSamples ; r++) {
          of = rhoar[r] - xar[m] * costheta;
          psi = of - yar[n] * sintheta;
          if ((2 * fabs(psi * costheta) < DeltaX) && (2 * fabs(psi * sintheta) < DeltaX))
            tempV[t + r * ThetaSamples] = DeltaX;
          else
            tempV[t + r * ThetaSamples] = 0.0;
        }
      }
      break;

    case _RNN: /* ray driven nearest neighbor */
      for (t = 0 ; t < ThetaSamples ; t++) {
        sintheta = sinthetaar[t];
        costheta = costhetaar[t];
        for (r = 0 ; r < RhoSamples ; r++) {
          if (sintheta > onedivsq2) {
            realn = (-Ymin + (rhoar[r] - Xmin * costheta) / sintheta) / DeltaY;
            intn = (int) floor(realn + 0.5);
            if ((intn >= 0) && (intn < YSamples))
              tempV[t + r * ThetaSamples] = DeltaX;
            realn -= costheta / sintheta;
          } else {
            realm = (-Xmin + (rhoar[r] - Ymin * sintheta) / costheta) / DeltaX;
            intm = (int) floor(realm + 0.5);
            if ((intm >= 0) && (intm < XSamples))
              tempV[t + r * ThetaSamples] = DeltaX;
          }
        }
      }
      break;

    case _RL: /* ray driven linear interpolation */
      for (t = 0 ; t < ThetaSamples ; t++) {
        sintheta = sinthetaar[t];
        costheta = costhetaar[t];
        if (sintheta > onedivsq2) {
          cst1 = Ymin + Xmin * costheta;
          cst2 = 1.0 / (sintheta * DeltaY);
        } else {
          cst1 = Xmin + Ymin * sintheta;
          cst2 = 1.0 / (costheta * DeltaX);
        }
        for (r = 0 ; r < RhoSamples ; r++) {
          if (sintheta > onedivsq2) {
            realn = (rhoar[r] - cst1) * cst2;
            intn = (int) floor(realn);
            if ((intn >= 0) && (intn < YSamples - 1)) {
              tempV[t + r * ThetaSamples] = DeltaX * (1 - (deltan = realn - intn));
              tempV[t + (r + 1) * ThetaSamples] = DeltaX * deltan;
            }
          } else {
            realm = (rhoar[r] - cst1) * cst2;
            intm = (int) floor(realm);
            if ((intm >= 0) && (intm < XSamples - 1)) {
              tempV[t + r * ThetaSamples] = DeltaX * (1 - (deltam = realm - intm));
              tempV[(t + 1) + r * ThetaSamples] = DeltaX * deltam;
            }
          }
        }
      }
      break;

    case _P1:        /* Uses same kernel */
    case _P2:        /* pixel driven interpolation for all points */
      for (t = 0 ; t < ThetaSamples ; t++)
        for (r = 0 ; r < RhoSamples ; r++)
          tempV[t + r * ThetaSamples] = SquareRadon(rhoar[r], thetaar[t], xar[m], yar[n]);
      break;
    default :Error("Wrong interpolation type detected");
  }
  return NewVector;
}

/***********************************************************************
[NAME]
MaskVec

[SYNOPSIS]
void MaskVec(Vector *MyVec,Vector *MyMask,
             float LowLimit,float HighLimit);

[DESCRIPTION] 

This function is used to mask an image {\tt MyImage} with another
image {\tt MyMask}. This can be used between in the reconstruction for
regularization. The two Vectors much have identical length. If the
sample in {\tt MyMask} is -1, and the pixel in {\tt MyVec} is set to
{\tt LowLimit}, if 1 to {\tt HighLimit}, and 0 implies that the sample
value is unchanged.

[USAGE]
{\tt MaskVec(MyVec,MyMask,0,10);}

Masks {\tt MyVec} with {\tt MyMask} and uses 0 for the low level, and 10
for the high level.

[REVISION]
March 29 1996 PT
***********************************************************************/
void MaskVec(Vector *MyVec, Vector *MyMask,
             float LowLimit, float HighLimit) {
  int n, N;
  float *ima, *mask;

  Print(_DDebug, "MaskVec\n");

  if ((N = MyVec->N) != MyMask->N)
    Error("MaskVec: Incompatible Image sizes");

  ima = MyVec->value;
  mask = MyMask->value;
  for (n = 0 ; n < N ; n++) {
    if (mask[n] == -1)
      ima[n] = LowLimit;
    else if (mask[n] == 1)
      ima[n] = HighLimit;
    else if (mask[n] != 0)
      Error("MaskVec: Invalid value detected in Mask");
  }
}
