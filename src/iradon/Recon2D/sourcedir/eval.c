/***************************************************************************

[HEADER] 

This file contains tools for evaluating images,
comparing different images them and calculating
some useful statistics about them.

***************************************************************************/

#include <math.h>
#include <stdio.h>
#include "imgtools.h"
#include "misc.h"
#include "calc.h"

/***************************************************************************
[NAME]
LargestValue

[SYNOPSIS]
float LargestValue(Image *MyImage)

[DESCRIPTION]

Finds the largest value in the image, and returns
it (for real images only).

[USAGE]
\tc{t=LargestValue(Test);}

Returns the largest value of \tc{Test} in \tc{t}. 

[REVISION]
Oct. 94, JJJ
***************************************************************************/
float LargestValue(Image *MyImage) {
  int m, n;
  float LV, *Data;

  for (LV = 0, m = 0 ; m < MyImage->M ; m++) {
    Data = MyImage->Signal[m];
    for (n = 0 ; n < MyImage->N ; n++)
      if (Data[n] > LV) {
        LV = Data[n];
      }
  }
  return LV;
}

/***************************************************************************
[NAME]
SmallestValue

[SYNOPSIS]
float SmallestValue(Image *MyImage)

[DESCRIPTION]

Finds the smallest value in the image, and
returns it (for real images only).

[USAGE]
\tc{t=SmallestValue(Test);}

Returns the smallest value of \tc{Test} in \tc{t}. 

[REVISION]
Oct. 94, JJJ
***************************************************************************/
float SmallestValue(Image *MyImage) {
  int m, n;
  float SV, *Data;

  for (SV = 0, m = 0 ; m < MyImage->M ; m++) {
    Data = MyImage->Signal[m];
    for (n = 0 ; n < MyImage->N ; n++)
      if (Data[n] <= SV)
        SV = Data[n];
  }
  return SV;
}

/***************************************************************************
[NAME]
MeanValue

[SYNOPSIS]
float MeanValue(Image *MyImage)

[DESCRIPTION]
Finds the mean value of the image (for real image only).

[USAGE]
\tc{t=MeanValue(Test);}

Returns the mean value of \tc{Test} in \tc{t}. 

[REVISION]
Oct. 94, JJJ
***************************************************************************/
float MeanValue(Image *MyImage) {
  int m, n;
  float *Data;
  float MV;

  for (MV = 0, m = 0 ; m < MyImage->M ; m++) {
    Data = MyImage->Signal[m];
    for (n = 0 ; n < MyImage->N ; n++)
      MV += Data[n];
  }
  MV /= (float) (MyImage->N * MyImage->M);
  return MV;
}

/***************************************************************************
[NAME]
Deviation

[SYNOPSIS]
float Deviation(Image *MyImage)

[DESCRIPTION] 

Finds and returns the standard deviation of the
image. Central estimator used (for real images only)

[USAGE]
\tc{t=Deviation(Test);}

Returns the deviation of \tc{Test} in \tc{t}. 

[REVISION]
Oct. 94, JJJ
***************************************************************************/
float Deviation(Image *MyImage) {
  int m, n;
  float MV, DV, *Data;

  MV = MeanValue(MyImage);
  for (DV = 0, m = 0 ; m < MyImage->M ; m++) {
    Data = MyImage->Signal[m];
    for (n = 0 ; n < MyImage->N ; n++)
      DV += (MV - Data[n]) * (MV - Data[n]);
  }
  DV = sqrt((float) DV / (MyImage->N * MyImage->M));
  return DV;
}

/***************************************************************************
[NAME]
L1Norm

[SYNOPSIS]
float L1Norm(Image *OrgImage, 
             Image *TestImage)

[DESCRIPTION]
Calculates the $L_1$--norm betveen the Original image and the test image.

[USAGE]
\tc{L1=L1Norm(RefImage,TestImage);}

returns the $L_1$ norm betveen {\tt RefImage} and {\tt TestImage}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
float L1Norm(Image *OrgImage, Image *TestImage) {
  int m, n;
  float Result, L1n, L1d, *data1, *data2, DC;

  if (!(OrgImage->M == TestImage->M) && (OrgImage->N == TestImage->N)
      && (OrgImage->ArrayType == TestImage->ArrayType))
    Error("Pictures must be same size (L1Norm)");

  L1n = 0.0;
  L1d = 0.0;
  DC = MeanValue(TestImage) - MeanValue(OrgImage);

  for (m = 0 ; m < OrgImage->M ; m++) {
    data1 = OrgImage->Signal[m];
    data2 = TestImage->Signal[m];
    for (n = 0 ; n < OrgImage->N ; n++) {
      L1n += (fabs(data2[n] - data1[n] + DC));
      L1d += (fabs(data1[n]));
    }
  }
  Result = L1n / L1d;
  return (Result);
}

/***************************************************************************
[NAME]
L2Norm

[SYNOPSIS]
float L2Norm(Image *OrgImage, 
             Image *TestImage)

[DESCRIPTION]
Calculates the $L_2$--norm betveen the Original image and the test image.

[USAGE]
\tc{L2=L2Norm(RefImage,TestImage);}

returns the $L_2$ norm betveen {\tt RefImage} and {\tt TestImage}.

[REVISION]
Oct. 94, JJJ and PAP\\
April 5, 96 PT bug - sign change on DC
***************************************************************************/
float L2Norm(Image *OrgImage, Image *TestImage) {
  int m, n;
  float L2, *data1, *data2, DC;

  if (!(OrgImage->M == TestImage->M) && (OrgImage->N == TestImage->N)
      && (OrgImage->ArrayType == TestImage->ArrayType))
    Error("Pictures must be same size (L2Norm)");

  L2 = 0.0;
  DC = MeanValue(TestImage) - MeanValue(OrgImage);

  for (m = 0 ; m < OrgImage->M ; m++) {
    data1 = OrgImage->Signal[m];
    data2 = TestImage->Signal[m];
    for (n = 0 ; n < OrgImage->N ; n++)
      L2 += ((data2[n] - data1[n] - DC) * (data2[n] - data1[n] - DC));
  }
  L2 = sqrt(1.0 / (OrgImage->M * OrgImage->N) * L2);
  L2 /= Deviation(OrgImage);
  return L2;
}

/***************************************************************************
[NAME]
DiffImage

[SYNOPSIS]
Image *DiffImage(Image *MyImage1,
                 Image *MyImage2);

[DESCRIPTION]

Returns an image containing the calculated difference between the two
input images calculated as \tc{MyImage1-MyImage2}. The images must be
equal in size.

[USAGE]
\tc{Diff=DiffImage(Image1,Image2);}

Returns the difference image in \tc{Diff}.

[REVISION]
Nov. 94, PAP
***************************************************************************/
Image *DiffImage(Image *MyImage1, Image *MyImage2) {
  Image *Diff;
  int m, n;
  float *Signal1, *Signal2, *DiffSignal;

  if ((MyImage1->M == MyImage2->M) && (MyImage1->N == MyImage2->N)
      && (MyImage1->ArrayType == MyImage2->ArrayType)) {
    Diff = NewFloatImage("DiffImage", MyImage1->M, MyImage1->N, MyImage1->ArrayType);
    for (m = 0 ; m < MyImage1->M ; m++) {
      Signal1 = MyImage1->Signal[m];
      Signal2 = MyImage2->Signal[m];
      DiffSignal = Diff->Signal[m];
      for (n = 0 ; n < MyImage1->N ; n++)
        DiffSignal[n] = Signal1[n] - Signal2[n];
    }
  } else {
    Diff = NULL;
    Error("Can not make subtraction of pictures\n");
  }
  return Diff;
}

/***************************************************************************
[NAME]
PrintStats

[SYNOPSIS]
void PrintStats(int DebugLevel,
                Image *MyImage)

[DESCRIPTION]

Prints some information about the image, some
statistics and the real world coordinates

[USAGE]
\tc{PrintStats(\_DNormal,Test);}

Prints statistics about \tc{Test}.

[REVISION]
Oct. 94, JJJ
***************************************************************************/
void PrintStats(int DebugLevel, Image *MyImage) {
  Print(DebugLevel,
        "Image stats for '%s', %s \n",
        MyImage->FileName,
        (MyImage->ArrayType == _RealArray) ? "real array" : "complex array");
  Print(DebugLevel,
        "  Parameters: Xmin = % 8.3e, DeltaX = %6.3e, NumX= %4d\n",
        MyImage->Xmin,
        MyImage->DeltaX,
        MyImage->M);
  Print(DebugLevel,
        "              Ymin = % 8.3e, DeltaY = %6.3e, NumY= %4d\n",
        MyImage->Ymin,
        MyImage->DeltaY,
        MyImage->N);
  Print(DebugLevel, "  Stats:    MinVal = % 6.2e, MaxVal = % 6.2e \n", SmallestValue(MyImage), LargestValue(MyImage));
  Print(DebugLevel, "            Mean   = % 6.2e, Dev.   = % 6.2e \n", MeanValue(MyImage), Deviation(MyImage));
}
