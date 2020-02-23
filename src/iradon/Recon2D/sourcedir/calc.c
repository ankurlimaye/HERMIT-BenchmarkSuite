/****************************************************************************** 
[HEADER]

The file calc.c contains the routines that does the actual
calculation, e.g. the transformation routines and the Fourier
transform routines.

11. Feb, 95 by JJJ and PAP.
12. Feb, 96 by PToft

******************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "imgtools.h"
#include "misc.h"
#include "calc.h"
#include "eval.h"

extern char *IniBuffer;

/******************************************************************************
[NAME]
Forward

[SYNOPSIS]
void Forward(Image *MyImage,Image *MySino)


[DESCRIPTION]

This function will forward project (Radon transform) the image {\tt MyImage}
into the sinogram {\tt MySino}.

[USAGE]
{\tt Forward(MyImage,MySino}

Forward projects the image {\tt MyImage} into the sinogram {\tt MySino}.


[REVISION]
Feb. 96, PToft
******************************************************************************/
void Forward(Image *MyImage, Image *MySino) {
  int t, r, m, n, mmin, mmax, nmin, nmax;
  int T, R, M, N;
  float sum, x_min, rhooffset, Delta_theta, Delta_rho, costheta, sintheta;
  float rho_min, theta, alpha, beta, nfloat, mfloat;
  float Delta_x, betap, eps;
  float **signal, **sino;

  eps = 1e-4;
  T = MySino->M;
  R = MySino->N;
  M = N = MyImage->M;

  rho_min = MySino->Ymin;
  Delta_theta = M_PI / T;
  Delta_rho = MySino->DeltaY;

  x_min = MyImage->Xmin;
  Delta_x = MyImage->DeltaX;
  signal = MyImage->Signal;
  sino = MySino->Signal;

  printf("Sino X=%i DX=%f Xmin=%f\n", MySino->M, MySino->DeltaX, MySino->Xmin);
  printf("Sino Y=%i DY=%f Ymin=%f\n", MySino->N, MySino->DeltaY, MySino->Ymin);
  printf("Image X=%i DX=%f Xmin=%f\n", MyImage->M, MyImage->DeltaX, MyImage->Xmin);
  printf("Image Y=%i DY=%f Ymin=%f\n", MyImage->N, MyImage->DeltaY, MyImage->Ymin);

  for (t = 0 ; t < T ; t++) {
    Print(_DNoLog, "Calculating: %4d of %d \r", t, T - 1);
    theta = t * Delta_theta;
    sintheta = sin(theta);
    costheta = cos(theta);
    rhooffset = x_min * (costheta + sintheta);
    if (sintheta > sqrt(0.5)) {
      alpha = -costheta / sintheta;
      for (r = 0 ; r < R ; r++) {
        beta = (r * Delta_rho + rho_min - rhooffset) / (Delta_x * sintheta);
        betap = beta + 0.5;
        sum = 0.0;
        if (alpha > 1e-6) {
          mmin = (int) ceil(-(beta + 0.5 - eps) / alpha);
          mmax = 1 + (int) floor((N - 0.5 - beta - eps) / alpha);
        } else if (alpha < -1e-6) {
          mmin = (int) ceil((N - 0.5 - beta - eps) / alpha);
          mmax = 1 + (int) floor(-(beta + 0.5 - eps) / alpha);
        } else {
          mmin = 0;
          mmax = M;
        }
        if (mmin < 0) mmin = 0;
        if (mmax > M) mmax = M;
        nfloat = betap + mmin * alpha;
        for (m = mmin ; m < mmax ; m++) {
          sum += signal[m][(int) nfloat];
          nfloat += alpha;
        }
        sino[t][r] = sum * Delta_x / fabs(sintheta);
      }
    } else {
      alpha = -sintheta / costheta;
      for (r = 0 ; r < R ; r++) {
        beta = (r * Delta_rho + rho_min - rhooffset) / (Delta_x * costheta);
        betap = beta + 0.5;
        sum = 0.0;
        if (alpha > 1e-6) {
          nmin = (int) ceil(-(beta + 0.5 - eps) / alpha);
          nmax = 1 + (int) floor((M - 0.5 - beta - eps) / alpha);
        } else if (alpha < -1e-6) {
          nmin = (int) ceil((M - 0.5 - beta - eps) / alpha);
          nmax = 1 + (int) floor(-(beta + 0.5 - eps) / alpha);
        } else {
          nmin = 0;
          nmax = N;
        }
        if (nmin < 0) nmin = 0;
        if (nmax > N) nmax = N;
        mfloat = betap + nmin * alpha;
        for (n = nmin ; n < nmax ; n++) {
          sum += signal[(int) mfloat][n];
          mfloat += alpha;
        }
        sino[t][r] = sum * Delta_x / fabs(costheta);
      }
    }
  }
  Print(_DNormal, "\n Finished         \n");
}

/******************************************************************************
[NAME]
BackProject

[SYNOPSIS]
void BackProject(Image *Sinogram, Image *InvMyImage)

[DESCRIPTION]

This function backprojects the sinogram {\tt Sinogram} onto the image
\tc{InvMyImage}. The sampling parameters in the inverse image
determines the actual projection geometry. This implementation uses
liniear interpolation in the $\rho$ direction to improve quality. An
{\tt Oversampling} Parameter can be specified in the {\tt INI}--file,
to further enhance the reconstruction quality.

[USAGE]
\tc{BackProject(TestSinogram,Test);}

Backprojects the sinogram \tc{TestSinogram} onto \tc{Test}.

[REVISION]
Oct. 94, JJJ and PAP\\
Feb. 95, JJJ, fixed minor bug.
******************************************************************************/
void BackProject(Image *Sinogram, Image *InvMyImage) {
  int m, n, N, t, XSamples, YSamples, RhoInt;
  Image *XCosTable, *YSinTable;
  float Xmin, Ymin, DeltaX, DeltaY, Rho, DeltaRho, RhoMax, TempSin, TempCos, sum;
  float *datam, *datan, *datat;

  Print(_DNormal, "Backprojecting sinogram...\n");
  PrintStats(_DDebug, Sinogram);
  PrintStats(_DDebug, InvMyImage);

  Xmin = InvMyImage->Xmin;
  Ymin = InvMyImage->Ymin;
  DeltaX = InvMyImage->DeltaX;
  DeltaY = InvMyImage->DeltaY;
  XSamples = InvMyImage->M;
  YSamples = InvMyImage->N;

  Print(_DDebug,
        "Transformed Image: min=(%.2f,%.2f), dx=%.2f dy=%.2f M=%d N=%d\n",
        Xmin, Ymin, DeltaX, DeltaY, XSamples, YSamples);

  /* Check if rho exceeded bounds in integration, if so stretch sinogram */
  RhoMax = sqrt(sq(max((fabs(Xmin)), (fabs(Xmin + (XSamples - 1) * DeltaX)))) +
      sq(max((fabs(Ymin)), (fabs(Ymin + (YSamples - 1) * DeltaY))))) + 1;
  if (RhoMax > (-Sinogram->Ymin)) {
    Print(_DDebug, "Required rho (%.2f) is bigger than maximum rho in sinogram (%.2f) \n", RhoMax, -Sinogram->Ymin);
    N = (int) (RhoMax / Sinogram->DeltaY) + 1;
    StretchImage(Sinogram, Sinogram->M, N * 2 + 1, _MiddleMiddle);
    PrintStats(_DDebug, Sinogram);
  }

  /* Allocate tabels and initialize */
  YSinTable = NewFloatImage("YSin", YSamples, Sinogram->M, _RealArray);
  XCosTable = NewFloatImage("XCos", XSamples, Sinogram->M, _RealArray);
  Print(_DDebug, "Initializing tables... \n");
  for (t = 0 ; t < Sinogram->M ; t++) {
    TempCos = cos(PI * t / Sinogram->M) / Sinogram->DeltaY;
    TempSin = sin(PI * t / Sinogram->M) / Sinogram->DeltaY;
    for (m = 0 ; m < XSamples ; m++)
      XCosTable->Signal[m][t] = ((float) m * DeltaX + Xmin) * TempCos;
    for (n = 0 ; n < YSamples ; n++)
      YSinTable->Signal[n][t] = ((float) n * DeltaY + Ymin) * TempSin - Sinogram->Ymin / Sinogram->DeltaY;
  }

  /* using linear interpolation to reduce noise in output image */
  Print(_DDebug, "Integrating... \n");
  for (m = 0 ; m < XSamples ; m++) {
    if (m % 5 == 0) Print(_DNoLog, "(%.2f pct. done) \r", (float) (m + 1) / XSamples * 100);
    datam = XCosTable->Signal[m];
    for (n = 0 ; n < YSamples ; n++) {
      datan = YSinTable->Signal[n];
      sum = 0.0;
      for (t = 0 ; t < Sinogram->M ; t++) {
        RhoInt = (int) (Rho = (datam[t] + datan[t]));
        /* if (Rho>=Sinogram->N || Rho<0) printf("* %d %d %d %f \n", m,n,t,Rho); */
        datat = &Sinogram->Signal[t][RhoInt];
        sum += datat[0] * (1 - (DeltaRho = (Rho - RhoInt)))
            + datat[1] * DeltaRho;
      }
      InvMyImage->Signal[m][n] = sum * Sinogram->DeltaX;
    }
  }
  FreeImage(XCosTable);
  FreeImage(YSinTable);
}

/******************************************************************************
[NAME]
IChirpSpectrum

[SYNOPSIS]
Image *IChirpSpectrum(Image *MyImage)

[DESCRIPTION]

This function calculates an inverse two dimensional Fourier transform
of the input image. The function uses the Chirp-z transformation
to zoom into the region of interest (ROI). The parametes of the
returned image are read from the \tc{.ini} file.

[USAGE]
\tc{Image=IChirpSpectrum(TestSpectrum);}

Inverse transforms \tc{TestSpectrum} into \tc{Image}.

[REVISION]
8. Nov. 94, JJJ
******************************************************************************/
Image *IChirpSpectrum(Image *MyImage) {
  int m;
  float rhomin, deltarho, rhosamples, phi0x, phi0y, theta0x, theta0y;
  Image *tempimage, *InvMyImage;

  Print(_DNormal, "Inverse transforming spectrum...\n");
  NormImage(MyImage, 1.0 / (IniFile.DeltaX * IniFile.DeltaY * MyImage->N * MyImage->M), 0.0);

  PrintStats(_DDebug, MyImage);
  deltarho = 1.0 / (2 * (float) MyImage->Xmin);
  rhosamples = (1 / MyImage->DeltaX) / deltarho;
  rhomin = -(rhosamples) / 2 * deltarho;
  Print(_DDebug, "RhoMin=%f numRho=%f, DeltaRho=%f\n", rhomin, rhosamples, deltarho);

  Print(_DDebug,
        "Transformed Image: (%3.2f,%3.2f), %3.2f %3.2f %d %d\n",
        IniFile.Xmin,
        IniFile.Ymin,
        IniFile.DeltaX,
        IniFile.DeltaY,
        IniFile.
            XSamples,
        IniFile.YSamples);

  phi0x = IniFile.DeltaX / ((rhosamples) * deltarho);
  phi0y = IniFile.DeltaY / ((rhosamples) * deltarho);
  theta0x = -0.5 * (IniFile.Xmin) / rhomin;
  theta0y = -0.5 * (IniFile.Ymin) / rhomin;
  Print(_DDebug, "px=%f py=%f tx=%f ty=%f \n", phi0x, phi0y, theta0x, theta0y);
  tempimage = NewFloatImage("temp", MyImage->N, IniFile.YSamples, _ComplexArray);

  /* chirp z inverse FFT */
  for (m = 0 ; m < MyImage->M ; m++) {
    Print(_DNoLog, "Chirp-Z (y): (%d:%d) \r", m + 1, MyImage->M);
    tempimage->Signal[m] = ComplexChirpZ(MyImage->Signal[m], MyImage->N, tempimage->N, theta0y, phi0y, _IFFT);
  }
  InvMyImage = NewFloatImage("invimage", IniFile.YSamples, IniFile.XSamples, _ComplexArray);
  MirrorImage(tempimage);

  for (m = 0 ; m < IniFile.YSamples ; m++) {
    Print(_DNoLog, "ChirpZ (x): (%d:%d)     \r", m + 1, tempimage->M);
    InvMyImage->Signal[m] = ComplexChirpZ(tempimage->Signal[m], tempimage->N, InvMyImage->N, theta0x, phi0x, _IFFT);
  }
  FreeImage(tempimage);
  MirrorImage(InvMyImage);
  RealImage(InvMyImage);
  InvMyImage->Xmin = IniFile.Xmin;
  InvMyImage->Ymin = IniFile.Ymin;
  InvMyImage->DeltaX = IniFile.DeltaX;
  InvMyImage->DeltaY = IniFile.DeltaY;
  PrintStats(_DDebug, InvMyImage);
  return (InvMyImage);
}

/******************************************************************************
[NAME]
IFFTSpectrum

[SYNOPSIS]
Image *IFFTSpectrum(Image *MyImage)

[DESCRIPTION]

This function calculates an inverse two dimensional Fourier transform
of the input image. The function uses the fast Fourier transformation
and cropping to zoom into the region of interest (ROI). The parametes
of the returned image are read from the \tc{.ini} file.

[USAGE]
\tc{Image=IFFTSpectrum(TestSpectrum);}

Inverse transforms \tc{TestSpectrum} into \tc{Image}.

[REVISION]
Jan. 95, JJJ and PAP
******************************************************************************/
Image *IFFTSpectrum(Image *MyImage) {
  int FFTLength;
  Image *InvMyImage;

  Print(_DNormal, "Inverse transforming spectrum...\n");
  PrintStats(_DDebug, MyImage);

  InvMyImage = CopyImage(MyImage);
  FFTLength = 1 << (int) (log(InvMyImage->M) / log(2) + 1);
  StretchImage(InvMyImage, FFTLength, FFTLength, _MiddleMiddle);
  FFTShift(InvMyImage, 2);

  FFTImage(InvMyImage, _IFFT);

  FFTShift(InvMyImage, 2);
  RealImage(InvMyImage);

  NormImage(InvMyImage, 1 / (IniFile.DeltaX * IniFile.DeltaY), 0.0);

  CropImage(InvMyImage, ((int) (IniFile.Xmin / IniFile.DeltaX + FFTLength / 2 + 0.5)),
            ((int) (IniFile.Ymin / IniFile.DeltaY + FFTLength / 2 + 0.5)), IniFile.XSamples, IniFile.YSamples);

  InvMyImage->Xmin = IniFile.Xmin;
  InvMyImage->Ymin = IniFile.Ymin;
  InvMyImage->DeltaX = IniFile.DeltaX;
  InvMyImage->DeltaY = IniFile.DeltaY;
  PrintStats(_DDebug, InvMyImage);
  return (InvMyImage);
}

/******************************************************************************
[NAME]
FilteredBack

[SYNOPSIS]
Image *FilteredBack(Image* MyImage)

[DESCRIPTION]
Calculates the discrete inverse Radon transformation of the sinogram
{\tt MyImage} using the filtered backprojection theorem. This function does the filtering, calls
\tc{BackProject} for the backprojection part.

[USAGE]
\tc{InvTest=FilteredBack(Test);}

Returns the inverse radon transformated image of \tc{Test} in \tc{InvTest}.

[REVISION]
Oct. 94, JJJ and PAP\\
Oct. 94, JJJ, moved the backprojection part into a seperat routine\\
Dec. 94, PAP, New filtering method.
******************************************************************************/
Image *FilteredBack(Image *MyImage) {
  int FFTLength, OldHeight;
  Image *InvMyImage;

  Print(_DNormal, "Filtered-Backproject transforming: '%s'\n", MyImage->FileName);
  Print(_DNormal, "Filtering...\n");

  Print(_DDebug, "Getting transformation parameters from INI file\n");

  /* Frequenzy filtering */
  OldHeight = MyImage->N;
  FFTLength = 1 << (int) (log(MyImage->N) / log(2) + 1 + 1);
  StretchImage(MyImage, MyImage->M, FFTLength, _LowerMiddle);
  RealVerticalFFT(MyImage, _FFT);
  FilterRealSpectrum(MyImage, IniFile.FilterType);
  StretchImage(MyImage, MyImage->M, MyImage->N * IniFile.InterPol, _LowerMiddle);
  MyImage->DeltaY /= IniFile.InterPol;
  RealVerticalFFT(MyImage, _IFFT);
  ShrinkImage(MyImage, MyImage->M, (OldHeight - 1) * IniFile.InterPol + 1, _LowerMiddle);

  /* Allocate new image, and put transformation parameters in it*/
  InvMyImage = NewFloatImage("InvMyImage", IniFile.XSamples, IniFile.YSamples, _RealArray);
  InvMyImage->Xmin = IniFile.Xmin;
  InvMyImage->Ymin = IniFile.Ymin;
  InvMyImage->DeltaX = IniFile.DeltaX;
  InvMyImage->DeltaY = IniFile.DeltaY;

  BackProject(MyImage, InvMyImage);
  PrintStats(_DDebug, InvMyImage);
  return InvMyImage;
}

/***************************************************************************
[NAME]
BackFilter

[SYNOPSIS]
Image *BackFilter(Image* MyImage)

[DESCRIPTION]

Calculates the discrete inverse Radon transformation of the sinogram
\tc{MyImage} using the filtering after backprojection theorem. Calls
\tc{BackProject} for the backprojection, and does the filtering
afterwards.

[USAGE]
\tc{InvTest=BackFilter(Test);}

Returns the inverse Radon transform of \tc{Test} in \tc{InvTest}.

[REVISION]
Oct. 94, PAP\\
Oct. 94, JJJ, Changed to use the {\tt BackProject} rutine\\
Dec. 94, PAP, New filtering method
***************************************************************************/
Image *BackFilter(Image *MyImage) {
  int i, m, n, mm, nn, M, N, OldHeight, OldWidth;
  int XSamples, YSamples, CenterM, CenterN;
  float Xmin, Ymin, Res, *Resx, *Resy, TempFloat, *TempPoint;
  Image *InvMyImage;

  Print(_DNormal, "Filter after Backproject transforming: '%s'\n", MyImage->FileName);
  Print(_DDebug, "Sinogram dimensions: M:%i N:%i\n", MyImage->M, MyImage->N);
  Print(_DNormal, "Filtering ...\n");

  M = N = (int) ((MyImage->N - 1) / (float) sqrt(2)) + 1;
  Print(_DDebug, "Backprojected image dim.: M:%i N:%i\n", IniFile.XSamples, IniFile.YSamples);

  OldHeight = IniFile.XSamples;
  OldWidth = IniFile.YSamples;
  XSamples = 1 << (int) (log(IniFile.XSamples) / log(2) + 1);
  YSamples = 1 << (int) (log(IniFile.YSamples) / log(2) + 1);
  Xmin = IniFile.Xmin + ((int) ((OldHeight - XSamples - 1) / 2)) * IniFile.DeltaX;
  Ymin = IniFile.Ymin + ((int) ((OldWidth - YSamples - 1) / 2)) * IniFile.DeltaY;

  /* Allocate new image, and put transformation parameters in it*/
  InvMyImage = NewFloatImage("InvMyImage", XSamples, YSamples, _RealArray);
  InvMyImage->Xmin = Xmin;
  InvMyImage->Ymin = Ymin;
  InvMyImage->DeltaX = IniFile.DeltaX;
  InvMyImage->DeltaY = IniFile.DeltaY;
  BackProject(MyImage, InvMyImage);

  /* Filter the backprojected image */

  FFTImage(InvMyImage, _FFT);
  CenterM = InvMyImage->M / 2;
  CenterN = InvMyImage->N / 2;

  Resx = FloatVector(InvMyImage->M);
  Resy = FloatVector(InvMyImage->N);

  for (m = 0 ; m < InvMyImage->M ; m++) {
    mm = m;
    if (mm >= CenterM) mm = InvMyImage->M - m;
    Resx[m] = (((float) mm) / InvMyImage->M) * (((float) mm) / InvMyImage->M) /
        (InvMyImage->DeltaX * InvMyImage->DeltaX);
    /*Print(_DDebug,"Resx %5.5f\n ",Resx[m]);*/
  }

  for (n = 0 ; n < InvMyImage->N ; n++) {
    nn = n;
    if (nn >= CenterN) nn = InvMyImage->N - n;
    Resy[n] = (((float) nn) / InvMyImage->N) * (((float) nn) / InvMyImage->N) /
        (InvMyImage->DeltaY * InvMyImage->DeltaY);
    /*Print(_DDebug,"Resy %5.5f\n",Resy[n]);*/
  }

  for (m = 0 ; m < InvMyImage->M ; m++) {
    TempFloat = Resx[m];
    TempPoint = InvMyImage->Signal[m];
    for (n = 0, i = 0 ; n < InvMyImage->N ; n++) {
      Res = sqrt(TempFloat + Resy[n]);
      /*Print(_DDebug,"Res %5.5f\n",Res);*/
      TempPoint[i++] *= Res;
      TempPoint[i++] *= Res;
    }
  }

  FFTImage(InvMyImage, _IFFT);
  ShrinkImage(InvMyImage, OldHeight, OldWidth, _MiddleMiddle);
  RealImage(InvMyImage);
  NormImage(InvMyImage, 1.0, -MeanValue(InvMyImage));
  PrintStats(_DDebug, InvMyImage);
  free(Resx);
  free(Resy);

  return (InvMyImage);
}

/******************************************************************************
[NAME]
CentralSliceNN

[SYNOSPSIS]
Image *CentralSliceNN(Image* MyImage)

[DESCRIPTION]

Calculates the discrete inverse Radon transformation of the sinogram
\tc{MyImage} using the central slice theorem with nearest neighbor
interpolation. Returns the interpolated spectrum.

[USAGE]
\tc{InvTest=CentralSliceNN(test);}

Returns the nearest neighbor interpolated spectrum of \tc{test} in \tc{InvTest}.

[REVISION]
Oct. 94, JJJ and PAP
******************************************************************************/
Image *CentralSliceNN(Image *MyImage) {
  int m, n, M, N, Realm, Realn, NyInt, ThetaInt, FFTLength, FFTLengthX, FFTLengthY, DistX, DistY;
  float DC, Ny, Theta, ThetaAdj, NySign, TempNy, TempTheta, NyMax, NyAdj;
  float *InvP, *MyInt0;
  Image *InvMyImage;

  Print(_DNormal, "Central Slice transforming '%s' \n", MyImage->FileName);

  if (IniFile.DeltaX != IniFile.DeltaY) Error("DeltaX must be equal to DeltaY\n");

  FFTLength = 1 << (int) (log(MyImage->N) / log(2) + 1);
  PrintStats(_DDebug, MyImage);

  /* Stretch, Shift and Fourier Transform image in rho direction */
  StretchImage(MyImage, MyImage->M, FFTLength, _MiddleMiddle);
  FFTShift(MyImage, 1);
  RealVerticalFFT(MyImage, _FFT);

  /* Expand image and place g(N/2) at the end of the vectors */
  /* And with one in theta direction to make interpolation easy!*/
  StretchImage(MyImage, MyImage->M + 1, MyImage->N + 1, _LowerLeft);
  for (m = 0 ; m < (MyImage->M - 1) ; m++) {
    MyImage->Signal[m][MyImage->N * 2 - 2] = MyImage->Signal[m][1];
    MyImage->Signal[m][1] = 0.0;
  }
  for (n = 0 ; n < MyImage->N ; n++) {
    MyImage->Signal[(MyImage->M - 1)][2 * n] = MyImage->Signal[0][2 * n];
    MyImage->Signal[(MyImage->M - 1)][2 * n + 1] = -MyImage->Signal[0][2 * n + 1];
  }

  DistX = (int) (2 * max(fabs(IniFile.Xmin), fabs(IniFile.Xmin + IniFile.XSamples * IniFile.DeltaX))
      / IniFile.DeltaX);
  DistY = (int) (2 * max(fabs(IniFile.Ymin), fabs(IniFile.Ymin + IniFile.YSamples * IniFile.DeltaY))
      / IniFile.DeltaY);
  FFTLengthX = 1 << (int) (log(DistX) / log(2) + 2);
  FFTLengthY = 1 << (int) (log(DistY) / log(2) + 2);

  /* New dimensions of image, must be odd  */
  M = FFTLengthX - 1;
  N = FFTLengthY - 1;
  InvMyImage = NewFloatImage("temp_spectrum", M, N, _ComplexArray);

  InvMyImage->Xmin = 0.5 / MyImage->DeltaY;
  InvMyImage->Ymin = 0.5 / MyImage->DeltaY;
  InvMyImage->DeltaX = 1.0 / (FFTLengthX * IniFile.DeltaX);
  InvMyImage->DeltaY = 1.0 / (FFTLengthY * IniFile.DeltaY);

  if ((FFTLengthX * IniFile.DeltaX < 1.42 * MyImage->DeltaY * 2 * (MyImage->N - 1)) ||
      (FFTLengthY * IniFile.DeltaY < 1.42 * MyImage->DeltaY * 2 * (MyImage->N - 1)))
    Print(_DNormal, "(Spectral information may be inadequate)\n");

  if ((IniFile.DeltaX > MyImage->DeltaY) || (IniFile.DeltaY > MyImage->DeltaY))
    Print(_DNormal, "(truncating spectrum)\n");

  Print(_DDebug, "New image:  M:%3d N:%3d \n", InvMyImage->M, InvMyImage->N);
  ThetaAdj = (MyImage->M - 1) / myPI;
  NyMax = 0.5 / MyImage->DeltaY;
  NyAdj = (MyImage->N - 1) / NyMax;

  /* Coordinate transform polar->rect */
  Print(_DNormal, "Transforming coordinates\n");
  for (m = 0, Realm = -(M - 1) / 2 ; m < M ; Realm++, m++) {
    InvP = InvMyImage->Signal[m];
    for (n = 0, Realn = -(N - 1) / 2 ; n < N ; Realn++, n++) {
      Ny = sqrt((float) Realn * (float) Realn * InvMyImage->DeltaX * InvMyImage->DeltaX +
          (float) Realm * (float) Realm * InvMyImage->DeltaY * InvMyImage->DeltaY);
      if (Ny < NyMax) {
        Theta = atan2((float) Realn, (float) Realm);

        if (Theta >= myPI) {
          Theta -= myPI;
          Ny *= -1;
        }
        if (Theta < 0) {
          Theta += myPI;
          Ny *= -1;
        }
        NySign = 1.0;
        if (Ny < 0) {
          Ny = -Ny;
          NySign *= -1;
        }

        TempNy = NyAdj * Ny;
        NyInt = (int) (TempNy + 0.5);

        TempTheta = ThetaAdj * Theta;
        ThetaInt = (int) (TempTheta + 0.5);

        MyInt0 = MyImage->Signal[ThetaInt];

        InvP[2 * n] = MyInt0[2 * NyInt];
        InvP[2 * n + 1] = MyInt0[2 * NyInt + 1];
        InvP[2 * n + 1] *= NySign;
      } else {
        InvP[2 * n] = 0.0;
        InvP[2 * n + 1] = 0.0;
      }
    }
  }

  for (m = 0, DC = 0.0 ; m < MyImage->M - 1 ; m++)
    DC += MyImage->Signal[m][0];
  InvMyImage->Signal[(InvMyImage->M - 1) / 2][(InvMyImage->N - 1)] = DC / (MyImage->M - 1);
  InvMyImage->Signal[(InvMyImage->M - 1) / 2][(InvMyImage->N - 1) + 1] = 0.0;
  NormImage(InvMyImage, MyImage->DeltaY, 0.0);
  PrintStats(_DDebug, InvMyImage);
  return (InvMyImage);
}

/******************************************************************************
[NAME]
CentralSliceBL

[SYNOSPSIS]
Image *CentralSliceBL(Image* MyImage)

[DESCRIPTION]

Calculates the discrete inverse Radon transformation of the sinogram
\tc{MyImage} using the central slice theorem with bi-linear
interpolation. Returns the interpolated spectrum.

[USAGE]
\tc{InvTest=CentralSliceBL(test);}

Returns the bi-linear interpolated spectrum of \tc{test} in \tc{InvTest}.

[REVISION] Oct. 94, JJJ and PAP
******************************************************************************/
Image *CentralSliceBL(Image *MyImage) {
  int
      m, n, M, N, Realm, Realn, NyInt, ThetaInt, FFTLength, FFTLengthX, FFTLengthY, DistX, DistY;
  float
      DC, Ny, Theta, ThetaAdj, NySign, DeltaNy, DeltaTheta, TempNy, TempTheta, NyMax, NyAdj;
  float *InvP, *MyInt0, *MyInt1;
  Image *InvMyImage;

  Print(_DNormal, "Central Slice transforming '%s' \n", MyImage->FileName);

  if (IniFile.DeltaX != IniFile.DeltaY) Error("DeltaX must be equal to DeltaY\n");

  FFTLength = 1 << (int) (log(MyImage->N) / log(2) + 1);
  PrintStats(_DDebug, MyImage);

  /* Stretch, Shift and Fourier Transform image in rho direction */
  StretchImage(MyImage, MyImage->M, FFTLength, _MiddleMiddle);
  FFTShift(MyImage, 1);
  RealVerticalFFT(MyImage, _FFT);

  /* Expand image and place g(N/2) at the end of the vectors */
  /* And with one in theta direction to make interpolation easy!*/
  StretchImage(MyImage, MyImage->M + 1, MyImage->N + 1, _LowerLeft);
  for (m = 0 ; m < (MyImage->M - 1) ; m++) {
    MyImage->Signal[m][MyImage->N * 2 - 2] = MyImage->Signal[m][1];
    MyImage->Signal[m][1] = 0.0;
  }
  for (n = 0 ; n < MyImage->N ; n++) {
    MyImage->Signal[(MyImage->M - 1)][2 * n] = MyImage->Signal[0][2 * n];
    MyImage->Signal[(MyImage->M - 1)][2 * n + 1] = -MyImage->Signal[0][2 * n + 1];
  }

  DistX = (int) (2 * max(fabs(IniFile.Xmin), fabs(IniFile.Xmin + IniFile.XSamples * IniFile.DeltaX))
      / IniFile.DeltaX);
  DistY = (int) (2 * max(fabs(IniFile.Ymin), fabs(IniFile.Ymin + IniFile.YSamples * IniFile.DeltaY))
      / IniFile.DeltaY);
  FFTLengthX = 1 << (int) (log(DistX) / log(2) + 2);
  FFTLengthY = 1 << (int) (log(DistY) / log(2) + 2);

  /* New dimensions of image, must be odd  */
  M = FFTLengthX - 1;
  N = FFTLengthY - 1;
  InvMyImage = NewFloatImage("temp_spectrum", M, N, _ComplexArray);

  InvMyImage->Xmin = 0.5 / MyImage->DeltaY;
  InvMyImage->Ymin = 0.5 / MyImage->DeltaY;
  InvMyImage->DeltaX = 1.0 / (FFTLengthX * IniFile.DeltaX);
  InvMyImage->DeltaY = 1.0 / (FFTLengthY * IniFile.DeltaY);

  if ((FFTLengthX * IniFile.DeltaX < 1.42 * MyImage->DeltaY * 2 * (MyImage->N - 1)) ||
      (FFTLengthY * IniFile.DeltaY < 1.42 * MyImage->DeltaY * 2 * (MyImage->N - 1)))
    Print(_DNormal, "(Spectral information may be inadequate)\n");

  if ((IniFile.DeltaX > MyImage->DeltaY) || (IniFile.DeltaY > MyImage->DeltaY))
    Print(_DNormal, "(truncating spectrum)\n");

  Print(_DDebug, "New image:  M:%3d N:%3d \n", InvMyImage->M, InvMyImage->N);
  ThetaAdj = (MyImage->M - 1) / myPI;
  NyMax = 0.5 / MyImage->DeltaY;
  NyAdj = (MyImage->N - 1) / NyMax;

  /* Coordinate transform polar->rect */
  Print(_DNormal, "Transforming coordinates\n");
  for (m = 0, Realm = -(M - 1) / 2 ; m < M ; Realm++, m++) {
    InvP = InvMyImage->Signal[m];
    for (n = 0, Realn = -(N - 1) / 2 ; n < N ; Realn++, n++) {
      Ny = sqrt((float) Realn * (float) Realn * InvMyImage->DeltaX * InvMyImage->DeltaX +
          (float) Realm * (float) Realm * InvMyImage->DeltaY * InvMyImage->DeltaY);
      if (Ny < NyMax) {
        Theta = atan2((float) Realn, (float) Realm);

        if (Theta >= myPI) {
          Theta -= myPI;
          Ny *= -1;
        }
        if (Theta < 0) {
          Theta += myPI;
          Ny *= -1;
        }
        NySign = 1.0;
        if (Ny < 0) {
          Ny = -Ny;
          NySign *= -1;
        }

        TempNy = NyAdj * Ny;
        NyInt = (int) (TempNy);
        DeltaNy = TempNy - NyInt;
        TempTheta = ThetaAdj * Theta;
        ThetaInt = (int) (TempTheta);
        DeltaTheta = TempTheta - ThetaInt;

        MyInt0 = MyImage->Signal[ThetaInt];
        MyInt1 = MyImage->Signal[ThetaInt + 1];

        InvP[2 * n] = ((MyInt0[2 * NyInt]) * (1.0 - DeltaNy) +
            (MyInt0[2 * NyInt + 2]) * (DeltaNy)) * (1.0 - DeltaTheta) +
            ((MyInt1[2 * NyInt]) * (1.0 - DeltaNy) +
                (MyInt1[2 * NyInt + 2]) * (DeltaNy)) * DeltaTheta;

        InvP[2 * n + 1] = ((MyInt0[2 * NyInt + 1]) * (1.0 - DeltaNy) +
            (MyInt0[2 * NyInt + 3]) * DeltaNy) * (1.0 - DeltaTheta) +
            ((MyInt1[2 * NyInt + 1]) * (1.0 - DeltaNy) +
                (MyInt1[2 * NyInt + 3]) * DeltaNy) * DeltaTheta;

        InvP[2 * n + 1] *= NySign;

      } else {
        InvP[2 * n] = 0.0;
        InvP[2 * n + 1] = 0.0;
      }
    }
  }

  for (m = 0, DC = 0.0 ; m < MyImage->M - 1 ; m++)
    DC += MyImage->Signal[m][0];
  InvMyImage->Signal[(InvMyImage->M - 1) / 2][(InvMyImage->N - 1)] = DC / (MyImage->M - 1);
  InvMyImage->Signal[(InvMyImage->M - 1) / 2][(InvMyImage->N - 1) + 1] = 0.0;
  NormImage(InvMyImage, MyImage->DeltaY, 0.0);
  PrintStats(_DDebug, InvMyImage);
  return (InvMyImage);
}

/******************************************************************************
[NAME]
CentralSliceCZ

[SYNOSPSIS]
Image *CentralSliceCZ(Image* MyImage)

[DESCRIPTION]

Calculates the discrete inverse Radon transformation of the sinogram
\tc{MyImage} using a variation of the central slice backprojection
theorem. By the use of nonlinear sampling of the Radon domain with the
use of the Chirp-z algorithm, the interpolation can be reduced to
simple one dimensional linear interpolation compared to normal CS.
The function returns the two dimensional spectrum of the reconstruced
image.

[USAGE]
\tc{InvTest=CentralSliceCZ(test);}

Returns the spectrum of the inverse Radon transform of
\tc{test} in \tc{InvTest}.

[REVISION]
Nov. 94, JJJ
******************************************************************************/

Image *CentralSliceCZ(Image *MyImage) {
  int m, n, M, N, thetaindex, thetaint, thetasamples, rhosamples;
  float DC, deltatheta, theta, thetaadj, *realtheta, rhomin, deltarho;
  Image *InvMyImage, *tempimage, *fourierimage;

  Print(_DNormal, "Chirp-Z Central Slice transforming '%s'\n", MyImage->FileName);

  Real2ComplexImage(MyImage);
  /* Interpolating sinogram in theta */
  M = (MyImage->M % 4 != 0) ? (MyImage->M - (4 - (MyImage->M % 4))) : (MyImage->M);
  tempimage = NewFloatImage("tempimage", M + 1, MyImage->N, _ComplexArray);
  realtheta = FloatVector(M + 1);

  Print(_DNormal, "Rebinning sinogram...\n");
  for (m = 0 ; m < M ; m++) {
    if (m >= M / 4 * 3)
      theta = PI - atan(4 / (float) M * ((float) M - (float) m));
    else if (m > M / 2)
      theta = PI / 2 + atan(4 * ((float) m - M / 2) / (float) M);
    else if (m <= M / 4)
      theta = atan(4 / (float) M * (float) m);
    else if (m < M / 2)
      theta = PI / 2 - atan(4 * (M / 2 - (float) m) / (float) M);
    else theta = PI / 2;
    realtheta[m] = theta; /* save the real angles */
    thetaadj = theta / PI * MyImage->M;
    deltatheta = thetaadj - (thetaindex = (int) thetaadj);
    for (n = 0 ; n < MyImage->N ; n++) {
      tempimage->Signal[m][2 * n] = MyImage->Signal[thetaindex][2 * n] * (1.0 - deltatheta)
          + MyImage->Signal[thetaindex + 1][2 * n] * deltatheta;
      tempimage->Signal[m][2 * n + 1] = MyImage->Signal[thetaindex][2 * n + 1] * (1.0 - deltatheta)
          + MyImage->Signal[thetaindex + 1][2 * n + 1] * deltatheta;
    }
  }

  /* Using Chirp-Z to interpolate in rho */
  Print(_DNormal, "Using chirp-z to interpolate in rho...\n");
  InvMyImage = NewFloatImage("invmyimage", tempimage->M, tempimage->N, _ComplexArray);
  for (m = 0 ; m < tempimage->M ; m++) {
    Print(_DNoLog, "Chirp-Z: (%d:%d) \r", m, tempimage->M);
    if ((m <= MyImage->M / 4) || (m > MyImage->M / 4 * 3)) {
      InvMyImage->Signal[m] = ComplexChirpZ(tempimage->Signal[m], tempimage->N, tempimage->N,
                                            0.0, fabs(0.5 / (tempimage->N * cos(realtheta[m]))), _FFT);
      for (n = (int) abs(cos(realtheta[m]) * tempimage->N) ; n < InvMyImage->N ; n++) {
        InvMyImage->Signal[m][2 * n] = 0;
        InvMyImage->Signal[m][2 * n + 1] = 0;
      }
    } else {
      InvMyImage->Signal[m] = ComplexChirpZ(tempimage->Signal[m], tempimage->N, InvMyImage->N,
                                            0.0, 0.5 / (tempimage->N * sin(realtheta[m])), _FFT);
      for (n = (int) ((sin(realtheta[m])) * InvMyImage->N) ; n < InvMyImage->N ; n++) {
        InvMyImage->Signal[m][2 * n] = 0;
        InvMyImage->Signal[m][2 * n + 1] = 0;
      }
    }
  }
  /* Fill M'th place in the array with theta(0) to ease up the interpolation */
  for (n = 0 ; n < InvMyImage->N ; n++) {
    InvMyImage->Signal[InvMyImage->M - 1][n * 2] = InvMyImage->Signal[0][2 * n];
    InvMyImage->Signal[InvMyImage->M - 1][n * 2 + 1] = InvMyImage->Signal[0][2 * n + 1];
  }
  FreeImage(tempimage);

  /* New image matches number of rho's to make better interpolation */

  fourierimage = NewFloatImage("fourierimage", InvMyImage->N * 2.0 - 1, InvMyImage->N * 2.0 - 1, _ComplexArray);
  M = fourierimage->M;
  N = fourierimage->N;
  thetasamples = InvMyImage->M - 1;

  /* Change samplingparameters for the spectrum */
  deltarho = MyImage->DeltaY;
  rhosamples = InvMyImage->N * 2.0;
  rhomin = -(fourierimage->M - 1) / 2 * MyImage->DeltaY;
  fourierimage->Xmin = 0.5 / (float) deltarho;
  fourierimage->Ymin = 0.5 / (float) deltarho;
  fourierimage->DeltaX = 1.0 / (rhosamples * deltarho);
  fourierimage->DeltaY = 1.0 / (rhosamples * deltarho);

  /* convert to cartesian grid */
  Print(_DNormal, "Converting polar->rectangular...\n");

  /* 0<=theta<PI/4 */
  for (m = 1 ; m <= (M - 1) / 2 ; m++)
    for (n = 0 ; n < m ; n++) {
      theta = (float) n / (float) m * (thetasamples / 4);
      thetaint = (int) (floor(theta));
      deltatheta = theta - thetaint;
      fourierimage->Signal[m + (M - 1) / 2][2 * (n + (N - 1) / 2)] =
          InvMyImage->Signal[thetaint][2 * m] * (1 - deltatheta)
              + InvMyImage->Signal[thetaint + 1][2 * m] * (deltatheta);
      fourierimage->Signal[m + (M - 1) / 2][2 * (n + (N - 1) / 2) + 1] =
          InvMyImage->Signal[thetaint][2 * m + 1] * (1 - deltatheta)
              + InvMyImage->Signal[thetaint + 1][2 * m + 1] * (deltatheta);
    }
  /* PI/4<=theta<3PI/4 */
  for (n = 1 ; n <= (N - 1) / 2 ; n++)
    for (m = -n + 1 ; m <= n ; m++) {
      theta = -(float) m / (float) n * (thetasamples / 4);
      thetaint = (int) (floor(theta));
      deltatheta = theta - thetaint;
      thetaint += +thetasamples / 2;
      fourierimage->Signal[m + (M - 1) / 2][2 * (n + (N - 1) / 2)] =
          InvMyImage->Signal[thetaint][2 * n] * (1 - deltatheta)
              + InvMyImage->Signal[thetaint + 1][2 * n] * (deltatheta);
      fourierimage->Signal[m + (M - 1) / 2][2 * (n + (N - 1) / 2) + 1] =
          InvMyImage->Signal[thetaint][2 * n + 1] * (1 - deltatheta)
              + InvMyImage->Signal[thetaint + 1][2 * n + 1] * (deltatheta);
    }
  /* 3PI/4<=theta<PI */
  for (m = -(M - 1) / 2 ; m <= -1 ; m++)
    for (n = 1 ; n <= -m ; n++) {
      theta = (float) n / (float) m * (thetasamples / 4);
      thetaint = (int) (floor(theta));
      deltatheta = theta - thetaint;
      thetaint += thetasamples;
      fourierimage->Signal[m + (M - 1) / 2][2 * (n + (N - 1) / 2)] =
          InvMyImage->Signal[thetaint][-2 * m] * (1 - deltatheta)
              + InvMyImage->Signal[thetaint + 1][-2 * m] * (deltatheta);
      fourierimage->Signal[m + (M - 1) / 2][2 * (n + (N - 1) / 2) + 1] =
          InvMyImage->Signal[thetaint][-2 * m + 1] * (1 - deltatheta)
              + InvMyImage->Signal[thetaint + 1][-2 * m + 1] * (deltatheta);
    }

  /* The rest of the spectrum is A complex conjugated copy of the first */
  for (m = 0 ; m < fourierimage->M ; m++)
    for (n = 0 ; n <= (fourierimage->N - 1) / 2 ; n++) {
      if ((n < (N - 1) / 2) || (n == (N - 1) / 2 && m < (M - 1) / 2)) {
        fourierimage->Signal[m][n * 2] =
            fourierimage->Signal[fourierimage->M - m - 1][(fourierimage->N - n - 1) * 2];
        fourierimage->Signal[m][n * 2 + 1] =
            -fourierimage->Signal[fourierimage->M - m - 1][(fourierimage->N - n - 1) * 2 + 1];
      }
    }

  /* Calc DC as the mean value for all theta */
  for (m = 0, DC = 0 ; m < thetasamples ; m++)
    DC += InvMyImage->Signal[m][0];
  fourierimage->Signal[(M - 1) / 2][(N - 1)] = DC / thetasamples;
  NormImage(fourierimage, IniFile.DeltaX * IniFile.DeltaY / MyImage->DeltaY, 0.0);

  return fourierimage;
}

/******************************************************************************
[NAME]
FilterRealSpectrum

[SYNOPSIS]
void FilterRealSpectrum(Image *MyImage,
                        int FilterType);

[DESCRIPTION]

Multiplies the image in the $y$ direction with a filter specified by
\tc{FilterType}. The image is treated like a spectrum from
\tc{RealFFT}. with only the frequencies from 0 to $f_s/2$ present. The
filters are the applied with the filter centrum on the first entry in
the image. The purpose of this routine is to do filtering before
backprojection. Filters are available:

\begin{tabbing}
\tc{FilterType} \== \tc{\_InvTriangle}\\
\>= \tc{\_Triangle}\\
\>= \tc{\_Hamming}\\
\>= \tc{\_Hanning}
\end{tabbing}

[USAGE]
\tc{FilterRealSpectrum(Test,\_Hamming);}

Multiplies the image \tc{Test} in the $y$ direction with a Hamming window.

[REVISION]
Oct. 94, JJJ\\
April 22 96 PT Cleanup
******************************************************************************/
void FilterRealSpectrum(Image *MyImage, int FilterType) {
  int m, n, N;
  float *Filter, *data;;

  Print(_DDebug, "FilterImage: Starting ... \n");

  if (MyImage->ArrayType != _ComplexArray)
    Error("Wrong Array format or FilterMode");
  N = MyImage->N + 1;

  Filter = FloatVector(N);

  for (n = 0 ; n <= MyImage->N ; n++) {
    switch (FilterType) {
      case _Ramp:
        if (n <= IniFile.FilterCutoff * N)
          Filter[n] = (float) n / (N * 2 * MyImage->DeltaY);
        else
          Filter[n] = 0.0;
        break;
      case _Hanning:
        if (n <= IniFile.FilterCutoff * N)
          Filter[n] = (float) n / (N * 2 * MyImage->DeltaY) *
              (0.5 + 0.5 * cos(IniFile.FilterCutoff * (float) n * myPI / (float) N));
        else
          Filter[n] = 0.0;
        break;
      case _Hamming:
        if (n <= IniFile.FilterCutoff * N)
          Filter[n] = (float) n / (N * 2 * MyImage->DeltaY) *
              (0.54 + 0.46 * cos(IniFile.FilterCutoff * (float) n * myPI / (float) N));
        else
          Filter[n] = 0.0;
        break;
      default:Error("FilterType unknown");
    }
  }
  for (m = 0 ; m < MyImage->M ; m++) {
    data = MyImage->Signal[m];
    for (n = 1 ; n < MyImage->N ; n++) {
      data[n * 2] *= Filter[n];
      data[n * 2 + 1] *= Filter[n];
    }
    data[0] *= Filter[0];
    data[1] *= Filter[MyImage->N];
  }
}

/******************************************************************************
[NAME]
FilterRealSpatial

[SYNOPSIS]
void FilterRealSpatial(Image *MyImage,
                       int FilterType);

[DESCRIPTION]

Multiplies the image in the $y$ direction with a filter specified by
\tc{FilterType}. The image is treated like a time sequence, e.g. the
filter centrum is placed in $N/2$. The purpose of this routine is to
do filtering (windowing) before Fourier transformation. These filters
are available

\begin{tabbing}
\tc{FilterType} \== \tc{\_InvTriangle}\\
\>= \tc{\_Triangle}\\
\>= \tc{\_Hamming}\\
\>= \tc{\_Hanning}
\end{tabbing}

[USAGE]
\tc{FilterRealSpatial(Test,\_Hanning);}

Multiplies the image \tc{Test} in the $y$ direction with a Hanning window.

[REVISION]
Oct. 94, JJJ\\
April 23 PT Removed triangle
******************************************************************************/
void FilterRealSpatial(Image *MyImage, int FilterType) {
  int m, n, N;
  float *Filter, *data;;

  Print(_DDebug, "FilterRealSpatial: Starting ... \n");

  if (MyImage->ArrayType != _RealArray)
    Error("Wrong Array format or FilterMode");
  N = MyImage->N;

  Filter = FloatVector(N);

  for (n = 0 ; n < N ; n++) {
    switch (FilterType) {
      case _Ramp:
        if (n >= N / 2) Filter[n] = 2 * ((float) n - N / 2) / N;
        else Filter[n] = 2 * (N / 2 - (float) n) / N;
        break;
      case _Hanning:Filter[n] = 0.5 * (1 - cos(2 * (float) n * myPI / (N)));
        break;
      case _Hamming:Filter[n] = 0.54 - 0.46 * cos(2 * (float) n * myPI / (N));
        break;
      default:Error("FilterType unknown");
    }
  }
  for (m = 0 ; m < MyImage->M ; m++) {
    data = MyImage->Signal[m];
    for (n = 0 ; n < MyImage->N ; n++)
      data[n] *= Filter[n];
  }
  Free(Filter);
}

/********************************************************************
[NAME]
ComplexChirpZ

[SYNOPSIS]
float* ComplexChirpZ(float* data,
                     int N,
                     int M,
                     float theta0,
                     float phi0,
                     int isign)

[DESCRIPTION]
This function will calculate the discrete spectrum in
frequency points linearly spaced

\begin{equation}
G(m)=\sum_{n=-\scbox{off}}^{N-1-\scbox{off}}
g(n)~e^{2\pi j\cdot \scbox{isign}\cdot n(m\cdot\phi_0+\theta_0)}
~m=0,1,...,M-1
\end{equation}

where $g(n)=\mbox{\tc{data}}[2*(n+\mbox{off})]~+~
~j~\mbox{\tc{data}}[2*(n+\mbox{off})+1]$
and $\mbox{off}=floor\left(\frac{N}{2}\right)$

The input vector inarr with length $2*\tc{N}$ contains the realpart
imaginary part of the first number (lowest \tc{n}), then real and
imaginary part of the second number (next \tc{n}), etc.  The vector
has the element corresponding to time 0 placed at the middle or just
right of the middle if \tc{N} is even.  The function returns an array
of length $2*\tc{M}$, with real and imaginary parts of the complex
spectrum. The spacing of the spectral components are controlled with
\tc{theta0}=$\theta_0$ and \tc{phi0}=$\phi_0>0$. If \tc{isign}=\_FFT
the spectrum is generated, and if \tc{isign}=\_IFFT the time signal is
generated from a given spectrum.

[USAGE]

\tc{arrout=ComplexChirpZ(inarr,40,80,0,1/80.0,\_FFT);}

This will generate a complex vector with 80 frequency points
spaced from 0 to half of the sampling frequency. The input vector
has 40 complex signal values.

[NOTE]

Remember to multiply the output vector with the spacing in frequency,
in order to get correct amplitude level.

[REVISION]
Nov. 94, JJJ and PT
********************************************************************/
float *ComplexChirpZ(float *data, int N, int M,
                     float theta0, float phi0, int isign) {
  float *Chirp, *har, *far;
  int l, L, m, n, off;
  float arg, piphi, twopitheta0;

  piphi = -PI * phi0 * isign;
  twopitheta0 = -2 * PI * theta0 * isign;

  L = (int) (1 << (int) (log(N + M - 1) / log(2.0) + 1));/* transformation length */

  Chirp = FloatVector(2 * M);
  far = FloatVector(2 * L);
  har = FloatVector(2 * L);

  off = (int) (N / 2);                           /* The zero time is in place 'off' */

  for (l = 0 ; l < N ; l++) {
    arg = -piphi * l * l - twopitheta0 * (l - off);
    far[0] = cos(arg);                        /* 'far' is used as temp array */
    far[1] = sin(arg);
    MultNew(far, &(data[2 * l]), &(har[2 * l]));  /* From 'data' to 'har' */
  }
  for (n = 0 ; n < M ; n++) {
    arg = piphi * n * n;
    far[2 * n] = cos(arg);
    far[2 * n + 1] = sin(arg);
  }
  for (n = L - N + 1 ; n < L ; n++) {
    arg = piphi * (n - L) * (n - L);
    far[2 * n] = cos(arg);
    far[2 * n + 1] = sin(arg);
  }
  for (m = 0 ; m < M ; m++)                         /* The return array is initialized */
  {
    arg = piphi * m * (2 * off - m);
    Chirp[2 * m] = cos(arg);
    Chirp[2 * m + 1] = sin(arg);
  }

  ComplexFFT(far, L, _FFT);                   /* FFT of far with length L */
  ComplexFFT(har, L, _FFT);                   /* FFT of har with length L */
  for (n = 0 ; n < L ; n++)
    MultReStore(&(far[2 * n]), &(har[2 * n]));   /* Multiply the two spectra -> far */
  ComplexFFT(far, L, _IFFT);                  /* IFFT of far gives far and har convoluted */
  for (m = 0 ; m < M ; m++)
    MultReStore(&(Chirp[2 * m]), &(far[2 * m])); /* Multiply with phase factors */
  free(har);
  free(far);
  return Chirp;                             /* Return Chirp */
}

/************************************************************************
[NAME]
ComplexZoomDFT

[SYNOSPSIS]
float* ComplexZoomDFT(float *data,
                      int nn,
                      int isign,
                      float theta0,
                      float phi0);

[DESCRIPTION]

This function calculates the chirp-z transformation (using the slow
DFT), of the vector \tc{data}. The function is equivalent to
\tc{ComplexChirpZ}, and may not be included in future releases.

[USAGE]
\tc{Test=ComplexZoomDFT(Test,100,100,0,1/100,\_IFFT);}

This does in fact calculate the inverse fourier transform of \tc{Test}.

[REVISION]
Oct. 94, JJJ and PT
**************************************************************************/
float *ComplexZoomDFT(float *data, int nn, float theta0, float phi0, int isign) {
  int n, m;
  float wpr, wpi, tempr;
  float *tempdata;

  tempdata = FloatVector(nn * 2);

  for (n = 0 ; n < nn ; n++) {
    for (m = 0 ; m < nn ; m++) {
      wpr = cos(2 * PI * m * (phi0 * n - theta0));
      wpi = sin(2 * PI * m * (phi0 * n - theta0));
      tempr = data[m * 2];
      tempdata[n * 2] += tempr * wpr - data[m * 2 + 1] * wpi;
      tempdata[n * 2 + 1] += tempr * wpi + data[m * 2 + 1] * wpr;
    }
  }
  Free(data);
  return tempdata;
}

/************************************************************************
[NAME]
ComplexFFT

[SYNOSPSIS]
void ComplexFFT(float *data,
                  int nn,
                  int isign);

[DESCRIPTION]

The function calculates the FFT of a complex sequence stored in
\tc{data}.  The array is an array of \tc{2*nn} float elements, i.e.,
\tc{nn} complex elements arranged like

data[0]=Re(g[0]), data[1]=Im(g[0])$,\ldots$\\
data[N-2]=Re(g[nn-1]), ar[N-1]=Im(g[nn-1]).

this sequence is transformed into the sequence

data[0]=Re(G[0]), ar[1]=Im(G[0])$,\ldots$\\
ar[N-2]=Re(G[nn-1]), ar[N-1]=Im(G[nn-1]).

where $N$ equals \tc{2*nn}. \tc{nn} must be a power of 2 (this is
checked). Isign can have the values \tc{\_FFT} or \tc{\_IFFT}.

This, and all other appearing fourier transform routines in this
library are modified versions of functions appearing in `Numerical
Recipes in C', 2'nd edition.

[USAGE]
 \tc{ComplexFFT(Test,128,\_FFT)}

Calculates the FFT of the sequence \tc{Test}, 128 complex numbers.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void ComplexFFT(float *data, int nn, int isign) {
  int n, mmax, m, j, istep, i;
  float wtemp, wr, wpr, wpi, wi, theta;
  float Temp, tempr, tempi;

/*  Print(_DDebug,"ComplexFFT: Starting ... \n");
*/
  n = nn << 1;
  j = 1;
  for (i = 1 ; i < n ; i += 2) {
    if (j > i) {
      Swap(data[j - 1], data[i - 1]);
      Swap(data[j], data[i]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  if (isign == _IFFT)
    for (i = 0 ; i < nn * 2 ; i++)
      data[i] /= nn;

  mmax = 2;
  while (n > mmax) {
    istep = 2 * mmax;
    theta = 6.28318530717959 / (isign * mmax);
    wpr = cos(theta);
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1 ; m < mmax ; m += 2) {
      for (i = m ; i <= n ; i += istep) {
        j = i + mmax;
        tempr = wr * data[j - 1] - wi * data[j];
        tempi = wr * data[j] + wi * data[j - 1];
        data[j - 1] = data[i - 1] - tempr;
        data[j] = data[i] - tempi;
        data[i - 1] += tempr;
        data[i] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi;
      wi = wi * wpr + wtemp * wpi;
    }
    mmax = istep;
  }
}

/************************************************************************
[NAME]
RealFFT

[SYNOSPSIS]
void RealFFT(float *data,
             int nn,
             int isign);

[DESCRIPTION]

The function calculates the FFT of a real sequence stored in
\tc{data}. The array is an array of \tc{nn} real elements. The
function uses \tc{ComplexFFT}. The data array for a time sequence must
be arranged like

data[0]=g[0], data[1]=g[1]$,\ldots,$data[N-2]=g[nn-1], ar[N-1]=g[nn-1],

where N equals \tc{nn}. A frequency sequence is returned as, and
must be arranged

data[0]=G[0], data[1]=G[N/2], data[2]=Re(G[1]), data[3]=Im(G[1])$,\ldots$\\
data[N-2]=Re(G[N/2-1]), data[N-1]=Im(G[N/2-1]).

Special care must be taken to ensure that data[1]=G[N/2] is in the right place.

\tc{nn} must be a power of 2. \tc{Isign} can have the values \tc{\_FFT} or
\tc{\_IFFT}.

[USAGE]
\tc{RealFFT(Test,128,\_FFT);}

Calculates the FFT of \tc{Test}, a real sequence of length 128.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void RealFFT(float *data, int n, int isign) {
  int i, i1, i2, i3, i4, np2;
  float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
  float wr, wi, wpr, wpi, wtemp, theta;

  Print(_DDebug, "RealFFT: Starting ... \n");

  theta = -myPI / (float) (n >> 1);

  if (isign == _FFT) {
    c2 = -0.5;
    ComplexFFT(data, n >> 1, _FFT);
  } else {
    c2 = 0.5;
    theta = -theta;
  }

  wpr = cos(theta);
  wpi = sin(theta);
  wr = wpr;
  wi = wpi;
  np2 = n + 2;
  for (i = 2 ; i <= (n >> 2) ; i++) {
    i4 = 1 + (i3 = np2 - 1 - (i2 = 1 + (i1 = i + i - 2)));
    h1r = c1 * (data[i1] + data[i3]);
    h1i = c1 * (data[i2] - data[i4]);
    h2r = -c2 * (data[i2] + data[i4]);
    h2i = c2 * (data[i1] - data[i3]);
    data[i1] = h1r + wr * h2r - wi * h2i;
    data[i2] = h1i + wr * h2i + wi * h2r;
    data[i3] = h1r - wr * h2r + wi * h2i;
    data[i4] = -h1i + wr * h2i + wi * h2r;
    wr = (wtemp = wr) * wpr - wi * wpi;
    wi = wi * wpr + wtemp * wpi;
  }
  if (isign == _FFT) {
    data[0] = (h1r = data[0]) + data[1];
    data[1] = h1r - data[1];
  } else {
    data[0] = c1 * ((h1r = data[0]) + data[1]);
    data[1] = c1 * (h1r - data[1]);
    ComplexFFT(data, n >> 1, _IFFT);
  }
}

void VerticalFFT(Image *MyImage, int isign) {
  Print(_DDebug, "VerticalFFT: Starting ... \n");
  if (MyImage->ArrayType == _RealArray)
    RealVerticalFFT(MyImage, isign);
  else
    ComplexVerticalFFT(MyImage, isign);
}

/************************************************************************
[NAME]
RealVerticalFFT

[SYNOSPSIS]
void RealVerticalFFT(image *MyImage,
                     int isign);

[DESCRIPTION]

This function calculates the FFT of all the vertical vectors in a real
image. If a FFT is to be calculated, is assumed that the image is
real. The image should be a complex spectrum if a IFFT is to be
preformed. The sequenzes should be stored like in \tc{RealFFT}. the
heigth of the image must be a power of two. This function is essential
a rewrite of \tc{RealFFT}, optimized for a lot of transformations of
the same length.

[USAGE]
\tc{RealFFT(Test,\_FFT);}

Calculates the FFT in the vertical direction for all of the image \tc{Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void RealVerticalFFT(Image *MyImage, int isign) {
  int i, n, mcount, i1, i2, i3, i4, np2;
  float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
  float wr, wi, wpr, wpi, wtemp, theta;
  float *data;

  if ((isign == _FFT) && (MyImage->ArrayType != _RealArray))
    Error("Image type (Complex/real) do not match transformation (FFT/IFFT)");

  if (isign == _FFT && MyImage->ArrayType == _RealArray) {
    n = MyImage->N;
    theta = -myPI / (float) (n >> 1);
    MyImage->ArrayType = _ComplexArray;
    MyImage->N /= 2;
    c2 = -0.5;
    ComplexVerticalFFT(MyImage, _FFT);
  } else {
    n = MyImage->N * 2;
    theta = myPI / (float) (n >> 1);
    c2 = 0.5;
  }

  wpr = cos(theta);
  wpi = sin(theta);
  wr = wpr;
  wi = wpi;
  np2 = n + 2;
  for (i = 2 ; i <= (n >> 2) ; i++) {
    i4 = 1 + (i3 = np2 - 1 - (i2 = 1 + (i1 = i + i - 2)));
    /* Step through Image one line at a time */
    for (mcount = 0 ; mcount < MyImage->M ; mcount++) {
      data = MyImage->Signal[mcount];
      h1r = c1 * (data[i1] + data[i3]);
      h1i = c1 * (data[i2] - data[i4]);
      h2r = -c2 * (data[i2] + data[i4]);
      h2i = c2 * (data[i1] - data[i3]);
      data[i1] = h1r + wr * h2r - wi * h2i;
      data[i2] = h1i + wr * h2i + wi * h2r;
      data[i3] = h1r - wr * h2r + wi * h2i;
      data[i4] = -h1i + wr * h2i + wi * h2r;
    }
    wr = (wtemp = wr) * wpr - wi * wpi;
    wi = wi * wpr + wtemp * wpi;
  }
  if (isign == _FFT) {
    for (mcount = 0 ; mcount < MyImage->M ; mcount++) {
      data = MyImage->Signal[mcount];
      data[0] = (h1r = data[0]) + data[1];
      data[1] = h1r - data[1];
    }
  } else {
    for (mcount = 0 ; mcount < MyImage->M ; mcount++) {
      data = MyImage->Signal[mcount];
      data[0] = c1 * ((h1r = data[0]) + data[1]);
      data[1] = c1 * (h1r - data[1]);
    }
    ComplexVerticalFFT(MyImage, _IFFT);
    MyImage->ArrayType = _RealArray;
    MyImage->N *= 2;
  }
}

/************************************************************************
[NAME]
ComplexVerticalFFT

[SYNOSPSIS]
void ComplexVerticalFFT(image *MyImage,
                          int isign);

[DESCRIPTION]

This function calculates the FFT of all the vertical vectors in a
complex image. The sequenzes must be stored like in
\tc{ComplexFFT}. the heigth of the image must be a power of two. This
function is essentially a rewrite of \tc{ComplexFFT}, optimized for a
lot of transformations of the same length.

[USAGE]
\tc{ComplexFFT(Test,\_FFT);}

Calculates the FFT in vertical direction for all of the image \tc{Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void ComplexVerticalFFT(Image *MyImage, int isign) {
  int n, nn, mmax, m, j, istep, i, mcount;
  float wtemp, wr, wpr, wpi, wi, theta;
  float Temp, tempr, tempi;
  float *data;

  if (MyImage->ArrayType == _RealArray)
    Error("Do not use 'ComplexVerticalFFT' for Real Images");

  CheckFFTLength(MyImage->N);
  nn = MyImage->N;
  n = nn << 1;
  j = 1;
  for (i = 1 ; i < n ; i += 2) {
    if (j > i) {
      /* Repeat for each line */
      for (mcount = 0 ; mcount < MyImage->M ; mcount++) {
        data = MyImage->Signal[mcount];
        Swap(data[j - 1], data[i - 1]);
        Swap(data[j], data[i]);
      }
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  if (isign == _IFFT)
    for (mcount = 0 ; mcount < MyImage->M ; mcount++) {
      data = MyImage->Signal[mcount];
      for (i = 0 ; i < nn * 2 ; i++)
        data[i] /= nn;
    }

  mmax = 2;
  while (n > mmax) {
    istep = 2 * mmax;
    theta = 6.28318530717959 / (isign * mmax);
    wpr = cos(theta);
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1 ; m < mmax ; m += 2) {
      /* Step through each line in Image one at a time */
      for (mcount = 0 ; mcount < MyImage->M ; mcount++) {
        data = MyImage->Signal[mcount];
        for (i = m ; i <= n ; i += istep) {
          j = i + mmax;
          tempr = wr * data[j - 1] - wi * data[j];
          tempi = wr * data[j] + wi * data[j - 1];
          data[j - 1] = data[i - 1] - tempr;
          data[j] = data[i] - tempi;
          data[i - 1] += tempr;
          data[i] += tempi;
        }
      }
      wr = (wtemp = wr) * wpr - wi * wpi;
      wi = wi * wpr + wtemp * wpi;
    }
    mmax = istep;
  }
}

/************************************************************************
[NAME]
ComplexNdimFFT

[SYNOSPSIS]
void ComplexVerticalFFT(float *data,
                        unsigned long *nn,
                        int ndim,
                        int isign);

[DESCRIPTION]

Calculates the n-dimensional Fourier transformation of a n-dimensional
data structure in data. (See Numerical Recipes pp. 523-524 for a
complete description). Used by \tc{FFTImage}. Data are stored in
sequence in \tc{data}. The dimensions in \tc{dim}, where \tc{dim[1]}
is the size in the first direction, and so on\ldots

[USAGE]
\tc{ComplexNdimFFT(Test,dim,3,\_FFT);}

Calculates the 3-dimensional FFT of the image stored in \tc{data}.

[REVISION]
Oct. 94, JJJ and PT
***************************************************************************/
void ComplexNdimFFT(float *data, unsigned long *nn, int ndim, int isign) {
  int idim, count;
  unsigned long i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
  unsigned long ibit, k1, k2, n, nprev, nrem, ntot;
  float tempi, tempr, Temp;
  double theta, wi, wpi, wpr, wr, wtemp;

  Print(_DDebug, "ComplexNDimFFT: Starting ...\n");
  for (ntot = 1, idim = 1 ; idim <= ndim ; idim++)
    ntot *= nn[idim];
  nprev = 1;
  for (idim = ndim ; idim >= 1 ; idim--) {
    n = nn[idim];
    nrem = ntot / (n * nprev);
    ip1 = nprev << 1;
    ip2 = ip1 * n;
    ip3 = ip2 * nrem;
    i2rev = 1;
    for (i2 = 1 ; i2 <= ip2 ; i2 += ip1) {
      if (i2 < i2rev) {
        for (i1 = i2 ; i1 <= i2 + ip1 - 2 ; i1 += 2) {
          for (i3 = i1 ; i3 <= ip3 ; i3 += ip2) {
            i3rev = i2rev + i3 - i2;
            Swap(data[i3], data[i3rev]);
            Swap(data[i3 + 1], data[i3rev + 1]);
          }
        }
      }
      ibit = ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
        i2rev -= ibit;
        ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1 = ip1;
    while (ifp1 < ip2) {
      ifp2 = ifp1 << 1;
      theta = isign * 6.28318530717959 / (ifp2 / ip1);
      wpr = cos(theta);
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (i3 = 1 ; i3 <= ifp1 ; i3 += ip1) {
        for (i1 = i3 ; i1 <= i3 + ip1 - 2 ; i1 += 2) {
          for (i2 = i1 ; i2 <= ip3 ; i2 += ifp2) {
            k1 = i2;
            k2 = k1 + ifp1;
            tempr = (float) wr * data[k2] - (float) wi * data[k2 + 1];
            tempi = (float) wr * data[k2 + 1] + (float) wi * data[k2];
            data[k2] = data[k1] - tempr;
            data[k2 + 1] = data[k1 + 1] - tempi;
            data[k1] += tempr;
            data[k1 + 1] += tempi;
          }
        }
        wr = (wtemp = wr) * wpr - wi * wpi;
        wi = wi * wpr + wtemp * wpi;
      }
      ifp1 = ifp2;
    }
    nprev *= n;
  }

  if (isign == _IFFT) {
    for (count = 1 ; count <= ntot * 2 ; count++) {
      data[count] /= ntot;
    }
  }
}

/************************************************************************
[NAME]
FFTImage

[SYNOSPSIS]
void FFTImage(Image *MyImage,
              int isign);

[DESCRIPTION]

Calculates the two dimensional discrete Fourier transformation of the image
using \tc{ComplexNdimFFT}. Appropriate checks are made to check if the
image can be transformed. Uses \tc{ComplexNdimFFT}.

[USAGE]
\tc{FFTImage(Test,\_IFFT);}

Calculates two dimensional FFT of the image stored in \tc{Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void FFTImage(Image *MyImage, int Isign) {
  float *Array;
  int m, n;
  unsigned long Dim[3];

  Print(_DDebug, "FFTImage: Starting ...\n");
  Array = FloatVector(2 * MyImage->M * MyImage->N);

  /* Convert Real image to Complex */
  if (MyImage->ArrayType == _RealArray) {
    StretchImage(MyImage, MyImage->M, MyImage->N * 2, _LowerMiddle);
    for (m = 0 ; m < MyImage->M ; m++)
      for (n = (MyImage->N / 2) - 1 ; n >= 0 ; n--) {
        MyImage->Signal[m][n * 2] = MyImage->Signal[m][n];
        MyImage->Signal[m][n * 2 + 1] = 0;
      }
    MyImage->N /= 2;
    MyImage->ArrayType = _ComplexArray;
  }

  /* Copy Image data form Image to array */
  for (m = 0 ; m < MyImage->M ; m++)
    memcpy(&Array[2 * m * MyImage->N], MyImage->Signal[m],
           2 * sizeof(float) * MyImage->N);

  Dim[1] = MyImage->M;
  Dim[2] = MyImage->N;
  CheckFFTLength(Dim[1]);
  CheckFFTLength(Dim[2]);

  Print(_DNormal, "FFTImage: FFT'ing '%s' in 2 dimensions\n", MyImage->FileName);

  ComplexNdimFFT(&Array[-1], Dim, 2, Isign);

  /* Copy Image data from array to Image */
  for (m = 0 ; m < MyImage->M ; m++)
    memcpy(MyImage->Signal[m], &Array[2 * m * MyImage->N],
           2 * sizeof(float) * MyImage->N);

  free(Array);
}

/************************************************************************
[NAME]
CheckFFTLength

[SYNOSPSIS]
void CheckFFTLength(int FFTLength)

[DESCRIPTION]

Checks if the \tc{FFTLength} is valid, i.e. a power of two. If not,
the program exits with an error message.

[USAGE]
\tc{CheckFFTLength(128);}

Checks if 128 is a valid FFTLength (it might be!).

[REVISION]
Oct. 94, JJJ
***************************************************************************/
void CheckFFTLength(int FFTLength) {
  int i;
  int isok;

  for (i = 2, isok = 0 ; (i <= 16384 && !isok) ; isok = (FFTLength == i), i *= 2);
  if (!isok)
    Error("FFTlength=%d is not a valid length\n", FFTLength);
  Print(_DDebug, "CheckFFTLength: Length of %d is ok.\n", FFTLength);
}

/************************************************************************
[NAME]
FFTShift

[SYNOSPSIS]
void FFTShift(Image *MyImage,
              int ShiftDimensions);

[DESCRIPTION]

Relocates the contents of the image. If \tc{ShiftDimensions}=1 the
image is shifted in the vertical direction so that the sequence

$-n,-n+1,\ldots-1,0,1,\ldots,n-1,n \Rightarrow  0,1,\ldots,n-1,n,-n,-n+1,\ldots,-1$

If \tc{ShiftDimensions}=2 the this is done is two directions. The
first and third quadrants are exchanged, as are the second and
fourth. The size of the image in the shifted dimension(s) must be an
even length.

[USAGE]
\tc{FFTShift(Test,2);}

Does a two dimensionel FFTShift on \tc{Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void FFTShift(Image *MyImage, int ShiftDimensions) {
  int m, Length, Middle;
  float *TempVector, *Temp;

  if (ShiftDimensions != 1 && ShiftDimensions != 2)
    Error("ShiftDimensions must be 1 or 2 (FFTShift)");

  if (MyImage->N % 2 || (MyImage->M % 2 && ShiftDimensions == 2))
    Error("Cannot FFTshift an odd array length");

  Print(_DDebug, "FFTShift: Shifting image '%s' in %d dimension(s)\n", MyImage->FileName, ShiftDimensions);

  TempVector = FloatVector(MyImage->N * MyImage->ArrayType);

  Length = MyImage->N * MyImage->ArrayType;
  Middle = (MyImage->N / 2) * MyImage->ArrayType;

  for (m = 0 ; m < MyImage->M ; m++) {
    memcpy(TempVector, &MyImage->Signal[m][Middle], (Length - Middle) * sizeof(float));
    memcpy(&TempVector[Length - Middle], MyImage->Signal[m], Middle * sizeof(float));
    Swap(TempVector, MyImage->Signal[m]);
  }

  if (ShiftDimensions == 2)
    for (m = 0 ; m < (MyImage->M) / 2 ; m++) {
      Swap(MyImage->Signal[m], MyImage->Signal[m + (MyImage->M) / 2]);
    }

  Free(TempVector);
}

void sort(int n, float *ra) {
  int l, j, ir, i;
  float rra;

  if (n < 2) return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1)
      rra = ra[--l - 1];
    else {
      rra = ra[ir - 1];
      ra[ir - 1] = ra[0];
      if (--ir == 1) {
        ra[0] = rra;
        return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && ra[j - 1] < ra[j]) ++j;
      if (rra < ra[j - 1]) {
        ra[i - 1] = ra[j - 1];
        j += (i = j);
      } else j = ir + 1;
    }
    ra[i - 1] = rra;
  }
}

/***************************************************************************
[NAME]
ImageFiltering

[SYNOPSIS]
ImageFiltering(Image *MyImage,
               int FilterType,
               unsigned int KernelX,
               unsigned int KernelY,
               unsigned int MeidianNumber)

[DESCRIPTION]

Filters the image with a filter of \tc{FilterType}. The two types of
filters are \\ \tc{\_Avarage}, calulate the average in the specified
kernel, does not use \tc{MedianNumber}. \\ \tc{\_Median}, replace with
the \tc{MedianNumber} counted from below.  \tc{MedianNumber} starts
from zero and have a maximum of \tc{KernelSize}.\\ \tc{KernelSize}
equals (2*KernelX+1)(2*KernelY+1), because KernelX specify number of
pixels on each side of the center pixel wich are inside the kernel in
x direction. KernelY is used as KernelX.

[REVISION]

Oct. 94, PAP\\
March 96 PT
***************************************************************************/
void ImageFiltering(Image *MyImage, int FilterType, int KernelX, int KernelY,
                    int MedianNumber) {
  int KernelPlace, Cerm, Cern, m, n, Count;
  int KernelWidthX, KernelWidthY, MaxNumber, nmin, nmax, mmin, mmax;
  float *KernelArea, Sum;
  Image *TempImage;

  if (MyImage->ArrayType == _ComplexArray)
    Error("Do not use a Complex Image in ImageFiltering\n");

  TempImage = CopyImage(MyImage);

  KernelWidthX = 2 * KernelX + 1;
  KernelWidthY = 2 * KernelY + 1;

  MaxNumber = KernelWidthX * KernelWidthY - 1;

  if (MedianNumber > MaxNumber)
    Print(_DNormal, "MediadNumber reduced to :%i\n", (MedianNumber = MaxNumber));

  if (!(KernelArea = (float *) malloc(KernelWidthX * KernelWidthY * sizeof(float))))
    Error("Memory allocation problems.\n");

  for (m = 0 ; m < MyImage->M ; m++) {
    mmin = max(m - KernelX, 0);
    mmax = min(m + KernelX + 1, MyImage->M - 1);
    for (n = 0 ; n < MyImage->N ; n++) {
      KernelPlace = 0.0;
      nmin = max(n - KernelY, 0);
      nmax = min(n + KernelY + 1, MyImage->N - 1);
      for (Cerm = mmin ; Cerm < mmax ; Cerm++)
        for (Cern = nmin ; Cern < nmax ; Cern++)
          KernelArea[KernelPlace++] = TempImage->Signal[Cerm][Cern];

      switch (FilterType) {
        case _Median: sort(KernelPlace, KernelArea);
          Count = (int) ((float) MedianNumber * ((float) KernelPlace / (float) MaxNumber));
          MyImage->Signal[m][n] = KernelArea[Count];
          break;
        case _Middle: Sum = 0.0;
          for (Count = 0 ; Count < KernelPlace ; Count++) Sum += KernelArea[Count];
          MyImage->Signal[m][n] = Sum / KernelPlace;
          break;
        default: Error("Wrong filter specified!!. \n");
      }
    }
  }
  free(KernelArea);
  FreeImage(TempImage);
}
