/***************************************************************************
[HEADER]

This file a group of two-dimensional image processing and manipulating
functions.

Some attention must be drawn to our {\tt .FIF} image format, it is a structure
defined as (direct copy from \tc{imgtools.h})

\begin{verbatim}
typedef struct {
  int  FIFIdType;        Id used to restore FIF:17737: _00EI           
  char FileName[100];    Name used for saving/restoring this image     
  char Description[80];  eg. short description of image (optional)     
  char Date[6];          Date (DDMMYY) (optional)                      
  int  Type;             Picture type (eg. format on disc, domain)     
  int  M;                Width of array                                
  int  N;                Height of array in no. numbers (Real/Complex) 
  int  ArrayType;        Defines number format: Complex(2) or Real(1)  
  float Xmin;            Leftmost coor. in original image              
  float Ymin;            Lowest coor. in original image                
  float DeltaX;          Quantisation distance in original image (X)   
  float DeltaY;          Quantisation distance in original image (Y)   
  float SignalMin;       Lowest signalvalue in array                   
  float SignalMax;       Highest signalvalue in array                  
  float ** Signal;       The raw data array                            
  } Image;								  
\end{verbatim}

if the image is real, the $(m,n)$'th value is stored in
\tc{Image->Signal[m][n]}. If the image is complex, the real part of the
complex number is stored in \tc{Image->Signal[m][2*n]}, and the
imaginary part in \tc{Image->Signal[m][2*n+1]}.

Last changed: 11. Nov, 94 by JJJ and PAP.

***************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "imgtools.h"
#include "gif_lib.h"
#include "analyze.h"
#include "misc.h"
#include "eval.h"

/***************************************************************************
[NAME]
NewFloatImage

[SYNOPSIS]
Image *NewFloatImage(char *FileName,
                     int M,
                     int N,
                     int ArrayType);

[DESCRIPTION] 

This function allocates and returns a float matrix of size {\tt M
$\times$ N}.  Initializes the variables {\tt N}, {\tt M}, {\tt
FileName} and \tc{Arraytype} in the Image structure.  {\tt ArrayType}
can be {\tt \_ComplexArray} or {\tt \_RealArray}.  The size of the
array is number of elements, either complex or real. The new image is
initialized to the value 0.

[USAGE] 

{\tt NewImage=NewFloatImage("test",7,8,\_RealArray);} 

Allocates a 7 $\times$ 8 matrix with real numbers.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
Image * NewFloatImage(char * FileName,int M,int N,int ArrayType)
{
  int m;
  Image * NewImage;

  Print(_DDebug,"NewFloatImage: Allocating memory for `%s'\n", FileName);

  if(!(NewImage=(Image *) malloc(sizeof(Image)))) 
    Error( "Memory allocation problems (NewFloatImage)" );

  if(!(NewImage->Signal=(float **) malloc(M*sizeof(float *)))) 
    Error("Memory allocation problems (NewFloatImage)" );

  for (m=0;m<M;m++)
    NewImage->Signal[m]=FloatVector((ArrayType*N));
  
  strcpy(NewImage->FileName,FileName);
  NewImage->M=M;
  NewImage->N=N;
  NewImage->Xmin=0;
  NewImage->Ymin=0;
  NewImage->DeltaX=1;
  NewImage->DeltaY=1;

  NewImage->ArrayType=ArrayType;

  NewImage->FIFIdType=_00EI;
  
  return NewImage;
}


/***************************************************************************
[NAME]
FreeImage

[SYNOPSIS]
void FreeImage(Image *MyImage);

[DESCRIPTION]
This function frees the memory allocated for {\tt MyImage}.
							 
[USAGE]
{\tt FreeImage(Test);}

Frees the memory allocated for \tc{test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/

void FreeImage(Image * MyImage)
{
  int m;

  Print(_DDebug,"FreeImage: Freeing `%s'\n", MyImage->FileName);

  for (m=0;m<MyImage->M;m++) 
    Free(MyImage->Signal[m]);
  Free(MyImage->Signal);
  Free(MyImage);
}


/***************************************************************************
[NAME]
ZeroImage

[SYNOPSIS]
void ZeroImage(Image *MyImage);

[DESCRIPTION]
This function initializes an entire image to the value 0.
							 
[USAGE]
{\tt ZeroImage(Test);}

Sets all values in {\tt Test} to 0.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void ZeroImage(Image * MyImage)
{
  int m,n;

  Print(_DDebug,"ZeroImage: Zeroing `%s'\n", MyImage->FileName);

    for (m=0;m<MyImage->M;m++)
      for (n=0;n<(MyImage->N*MyImage->ArrayType);n++)
        MyImage->Signal[m][n]=0.0;
}


/***************************************************************************
[NAME]
NormImage

[SYNOPSIS]
void NormImage(Image *MyImage,
               float Alpha,
               float Beta);

[DESCRIPTION] 

Normalization of images. This function calculates the following
function for every pixel ($x$) in the image: $x=\alpha(x+\beta)$.

[USAGE]

\tc{NormImage(Test,0.5,-0.1);}

Calculates $\mathrm{Test}(m,n)=0.5(\mathrm{Test}(m,n)-0.1)$ for all
 values of $m$ and $n$ for the image {\tt Test}.  

[REVISION] 
Oct. 94, JJJ
***************************************************************************/
void NormImage(Image * MyImage, float Alpha, float Beta)
{
  int m,n;
  float *data;

  Print(_DDebug,"MultiplyImage: Normating  `%s' with alpha=%f, beta=%f\n",
	MyImage->FileName,Alpha,Beta);

  if (MyImage->ArrayType==_RealArray) {
    for (m=0;m<MyImage->M;m++) {
      data=MyImage->Signal[m];
      for (n=0;n<MyImage->N;n++)
	data[n]=(data[n]+Beta)*Alpha;
    }
  }
  else
    for (m=0;m<MyImage->M;m++) {
      data=MyImage->Signal[m];
      for (n=0;n<MyImage->N;n++) {
	data[n*2  ] = (data[2*n]+Beta)*Alpha;
	data[n*2+1]*= Alpha;
      }
    }
}


/***************************************************************************
[NAME]
RenameImage

[SYNOPSIS]
void RenameImage(Image *MyImage, 
                 char *NewFileName);

[DESCRIPTION]
Renames the variable {\tt MyImage->FileName} to {\tt NewFileName}.

[USAGE]
{\tt RenameImage(Test, "NewTest");}

Renames the image {\tt Test} to ``NewTest''.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void RenameImage(Image *MyImage, char *NewFileName)
{
  Print(_DDebug,"RenameImage: Renaming `%s' to `%s'\n",
          MyImage->FileName, NewFileName);
  strcpy(MyImage->FileName,NewFileName);
}


/***************************************************************************
[NAME]
InitImage

[SYNOPSIS]
void InitImage(Image *MyImage);

[DESCRIPTION]

this function initializes the `real world parameters' in the structure
{\tt Image}.  Can be used when reading pictures without this
information in them (e.g. GIF). The image is initialized like a
standard sinogram, with the $x$ coordinates ranging from 0 to $\pi$,
and the $y$ coordinates centered around 0, with $\Delta y=1$.  The
following variables are initialized: {\tt DeltaX, DeltaY, Xmin, Ymin,
SignalMin, SignalMax}.

[USAGE]
{\tt InitImage(Test);}

Initializes the sinogram {\tt Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void InitImage(Image * MyImage)
{
  Print(_DDebug,"InitImage: Initializing `%s'\n", MyImage->FileName);
  MyImage->Xmin      = 0;
  MyImage->DeltaX    = PI/MyImage->M;
  MyImage->DeltaY    = 1;
  MyImage->Ymin      = -(MyImage->N-1)/2;
  MyImage->SignalMin = (SmallestValue(MyImage));
  MyImage->SignalMax = (LargestValue(MyImage));
}


/***************************************************************************
[NAME]
CopyImage

[SYNOPSIS]
Image *CopyImage(Image *MyImage);

[DESCRIPTION]

Creates a new image and copies the contents of the old one into
it. All memory allocation is done internally. Returns a pointer to the
new image

[USAGE]
{\tt Test2 = CopyImage(Test);}

Copies the Image {\tt Test} into {\tt Test2}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
Image *CopyImage(Image * MyImage)
{
  int m;
  float **TempSignal;
  Image *NewImage;

  Print(_DDebug,"CopyImage: Copying `%s'\n", MyImage->FileName);

  NewImage=NewFloatImage(MyImage->FileName, MyImage->M, MyImage->N,
                         MyImage->ArrayType);

  /* Copy Image header */
  TempSignal=NewImage->Signal;
  memcpy(NewImage, MyImage, sizeof(Image));
  NewImage->Signal=TempSignal;

  for(m=0; m<MyImage->M; m++) 
    memcpy(NewImage->Signal[m],MyImage->Signal[m],
           MyImage->N*MyImage->ArrayType*sizeof(float));

  return(NewImage);
}


/***************************************************************************
[NAME]
MirrorImage

[SYNOPSIS]
void MirrorImage(Image *MyImage);

[DESCRIPTION]

Flips the contents of the image {\tt MyImage} around the axis
$m=n$. Useful for doing multi-dimensional FFT's with normal
one-dimensional FFT's.

[USAGE]
{\tt MirrorImage(Test);}

Mirrors the Image {\tt Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void MirrorImage(Image* MyImage)
{
  Image *TempImage;
  int m,n, Temp;
  float **TempSignal;
  float tempf;
  Print(_DDebug,"MirrorImage: Mirroring `%s' \n", MyImage->FileName);

  TempImage=NewFloatImage("TMirror",MyImage->N, MyImage->M, MyImage->ArrayType);

  /* Copy data int new array */
  if (MyImage->ArrayType==_RealArray) 
    for (m=0; m<MyImage->M; m++) 
      for (n=0; n<MyImage->N; n++) 
        TempImage->Signal[n][m]=MyImage->Signal[m][n];
  else
    for (m=0; m<MyImage->M; m++) 
      for (n=0; n<(MyImage->N); n++) {
        TempImage->Signal[n][2*m  ]=MyImage->Signal[m][2*n  ];
        TempImage->Signal[n][2*m+1]=MyImage->Signal[m][2*n+1];
      }

  /* Exchange Signals */
  TempSignal=TempImage->Signal;
  TempImage->Signal=MyImage->Signal;
  MyImage->Signal=TempSignal;

  /* Exchange coords. */
  Swap(MyImage->M,TempImage->M);
  Swap(MyImage->N,TempImage->N);
  tempf=MyImage->DeltaX;
  MyImage->DeltaX=MyImage->DeltaY;
  MyImage->DeltaY=tempf;
  tempf=MyImage->Xmin;
  MyImage->Xmin=MyImage->Ymin;
  MyImage->Ymin=tempf;
  
  FreeImage(TempImage);
}


/***************************************************************************
[NAME]
Real2ComplexImage

[SYNOPSIS]
void Real2ComplexImage(Image *MyImage);

[DESCRIPTION]

Converts (expands) an image from real to complex arrayformat. The
image size in units remains the same.

[USAGE]
{\tt Real2ComplexImage(Test);}

Converts \tc{Test} from real to complex.

[REVISION]
Nov. 94, JJJ.
***************************************************************************/
void Real2ComplexImage(Image *MyImage)
{  
  int m,n;
  /* Convert Real image to Complex */
  if (MyImage->ArrayType==_RealArray) {
    StretchImage(MyImage, MyImage->M, MyImage->N*2, _LowerMiddle);
    for(m=0; m<MyImage->M; m++)
      for(n=(MyImage->N/2)-1; n>=0; n--) {
        MyImage->Signal[m][n*2]   = MyImage->Signal[m][n];
        MyImage->Signal[m][n*2+1] = 0;
      }
    MyImage->N/=2;
    MyImage->ArrayType=_ComplexArray;
  }
  Print(_DDebug,"Exiting Real2Complex\n");
}


/***************************************************************************
[NAME]
AbsoluteImage

[SYNOPSIS]
void AbsoluteImage(Image *MyImage);

[DESCRIPTION]

Calculates the absolute values of all the complex numbers in an image,
and returns the result in the same image. The image is converted from
complex to real format.

[USAGE]
{\tt AbsoluteImage(Test);}

Returns the numeric value of all elements of {\tt Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void AbsoluteImage(Image *MyImage)
{
  int m,n;

  Print(_DDebug,"AbsoluteImage: Absolut'ing `%s'\n", MyImage->FileName);

  if (MyImage->ArrayType==_ComplexArray) {

    for (m=0; m<MyImage->M; m++) 
      for (n=0; n<MyImage->N; n++) 
        MyImage->Signal[m][n]=
          sqrt(MyImage->Signal[m][n*2  ]*MyImage->Signal[m][n*2  ]+
               MyImage->Signal[m][n*2+1]*MyImage->Signal[m][n*2+1]);

    MyImage->ArrayType=_RealArray;
    MyImage->N=MyImage->N*2;

    ShrinkImage(MyImage,MyImage->M,MyImage->N/2,_LowerMiddle);
  }  
  else {
    for (m=0; m<MyImage->M; m++) 
      for (n=0; n<MyImage->N; n++) 
        MyImage->Signal[m][n]=fabs(MyImage->Signal[m][n]);
  }
}


/***************************************************************************
[NAME]
RealImage

[SYNOPSIS]
void RealImage(Image *MyImage);

[DESCRIPTION]

Finds the real values of a all the complex numbers in an image, and
returns the result in the same image. The image type is converted from
complex to real.

[USAGE]
{\tt RealImage(Test);}

Returns the real value of all members of {\tt Test}.

[REVISION]
Nov. 94, PAP
***************************************************************************/
void RealImage(Image *MyImage)
{
  int m,n;

  Print(_DDebug,"RealImage: Real'ing `%s'\n", MyImage->FileName);

  if (MyImage->ArrayType==_ComplexArray) {

    for (m=0; m<MyImage->M; m++) 
      for (n=0; n<MyImage->N; n++) 
        MyImage->Signal[m][n]=MyImage->Signal[m][n*2];
    MyImage->ArrayType=_RealArray;
    MyImage->N=MyImage->N*2;

    ShrinkImage(MyImage,MyImage->M,MyImage->N/2,_LowerMiddle);
  }  
}


/***************************************************************************
[NAME]
ImagImage

[SYNOPSIS]
void ImagImage(Image *MyImage);

[DESCRIPTION]

Finds the imaginary values of a all the complex numbers in an image,
and returns the result in the same image. The image type is converted from
complex to real.

[USAGE]
{\tt ImagImage(Test);}

Returns the imaginary value of all elements of {\tt Test}.

[REVISION]
Nov. 94, PAP
***************************************************************************/
void ImagImage(Image *MyImage)
{
  int m,n;

  Print(_DDebug,"ImagImage: Imag'ing `%s'\n", MyImage->FileName);

  if (MyImage->ArrayType==_ComplexArray) {

    for (m=0; m<MyImage->M; m++) 
      for (n=0; n<MyImage->N; n++) 
        MyImage->Signal[m][n]=MyImage->Signal[m][n*2+1];
    MyImage->ArrayType=_RealArray;
    MyImage->N=MyImage->N*2;

    ShrinkImage(MyImage,MyImage->M,MyImage->N/2,_LowerMiddle);
  }  
}


/***************************************************************************
[NAME]
CropImage

[SYNOPSIS]
void CropImage(Image *MyImage,
               int StartM,
               int StartN,
               int NewWidth,
               int NewHeigth);

[DESCRIPTION]

This function takes out a sub-image of the original image at the
coordinates \tc{(StartM,StartN)} with size {\tt NewWidth $\times$
NewHeight}. Appropriate tests are carried out to ensure that the
sub-image doesn't cross the edge of the original image. The sub-image
is returned in the original image.  i.e.  the original image is
destroyed.

[USAGE] 

{\tt CropImageImage(Test,10,10,25,25);}

Takes out a section of the size 25 $\times$ 25 starting with lower
left coordinates $(10,10)$ of the image {\tt Test}.

[REVISION]
Oct. 94, JJJ
***************************************************************************/
void CropImage(Image *MyImage, int NewStartM, int NewStartN, 
              int NewWidth, int NewHeight)
{
  int m, Temp;
  Image *TempImage;
  float **TempSignal;
  
  if ((NewStartM+NewWidth)>MyImage->M || ((NewStartN+NewHeight))>MyImage->N)
    Error("New Image is out of bounds (CropImage)");
 
  Print(_DDebug,"Cropping `%s' (%dx%d) to new dimensions (%dx%d)\n",
          MyImage->FileName, MyImage->M, MyImage->N, NewWidth, NewHeight);

  TempImage=NewFloatImage(MyImage->FileName,NewWidth,NewHeight,MyImage->ArrayType);

  /* update image real coords. */
  MyImage->Xmin=MyImage->Xmin+NewStartM*MyImage->DeltaX;
  MyImage->Ymin=MyImage->Ymin+NewStartM*MyImage->DeltaY;

  /* Copy Image Data */
  for(m=0;m<NewWidth;m++) {
    memcpy(TempImage->Signal[m],
           &MyImage->Signal[NewStartM+m][NewStartN*MyImage->ArrayType],
                            NewHeight*MyImage->ArrayType*sizeof(float));
  } 

  /* Exchange Signals */
  TempSignal=TempImage->Signal;
  TempImage->Signal=MyImage->Signal;
  MyImage->Signal=TempSignal;

  /* Exchange coords. */
  Swap(MyImage->M,TempImage->M);
  Swap(MyImage->N,TempImage->N);
  
  FreeImage(TempImage);
}


/***************************************************************************
[NAME]
StretchImage

[SYNOPSIS]
void StretchImage(Image *MyImage,
                  int NewWitdth,
                  int NewHeigth,
                  int NewArea);

[DESCRIPTION]

Stretches the image to new dimensions {\tt NewWidth $\times$
NewHeight}. The switch {\tt NewArea} determines where the old image
should be inserted in the new (and bigger) image. The new areas of the
image are initialized to 0.  {\tt NewArea} can have the following
values

\begin{verbatim}
NewArea = \_UpperLeft,  \_UpperMiddle,  \_UpperRight  
          \_MiddleLeft, \_MiddleMiddle, \_MiddleRight
          \_LowerLeft,  \_LowerMiddle,  \_LowerRight
\end{verbatim}

If case one of the `middle' options is used, the center of the image in
that (or both) directions will be placed according to the definition.

[USAGE]
{\tt StretchImage(Test, 102,102,\_MiddleMiddle);}

Stretches {\tt Test} to dimensions 102 $\times$ 102. The image is padded
in all directions e.g. the old image is placed in the middle.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void StretchImage(Image *MyImage, int NewWidth, int NewHeight, int NewArea)
{
  int NewStartM,NewStartN,m;
  Image *TempImage;
  int Temp;
  float **TempSignal;
 
  if(NewWidth<MyImage->M||NewHeight<MyImage->N)
    Error("Image size must be bigger than or equal to old size");

  TempImage=NewFloatImage(MyImage->FileName,
                          NewWidth,NewHeight,MyImage->ArrayType);

  switch(NewArea) {
  case _UpperLeft: case _MiddleLeft: case _LowerLeft:
   	NewStartM=0;	
      	break;
  case _UpperRight: case _MiddleRight: case _LowerRight:
	NewStartM=NewWidth-MyImage->M;
	break;
  case _UpperMiddle: case _MiddleMiddle: case _LowerMiddle:
	NewStartM=(NewWidth-MyImage->M)/2;
	if((MyImage->M%2)&&(!(NewWidth%2))) NewStartM++;
	break;
  default:
	NewStartM=0;
	Error("Unknown Area (StretchImage)");
  }
  switch(NewArea) {
  case _UpperLeft: case _UpperMiddle: case _UpperRight:
  	NewStartN=NewHeight-MyImage->N;	
      	break;
  case _MiddleLeft: case _MiddleMiddle: case _MiddleRight:
	NewStartN=(NewHeight-MyImage->N)/2;
	if((MyImage->N%2)&&(!(NewHeight%2))) NewStartN++;
	break;
  case _LowerLeft: case _LowerMiddle: case _LowerRight:
	NewStartN=0;
	break;
  default: 
        NewStartN=0;
	Error("Unknown Area (StretchImage)");
  }
  /* update image real coords. */
  MyImage->Xmin=MyImage->Xmin-NewStartM*MyImage->DeltaX;
  MyImage->Ymin=MyImage->Ymin-NewStartN*MyImage->DeltaY;
  
  Print(_DDebug,"StretchImage: Stretching `%s' (%dx%d) to new dimensions (%dx%d)\n",
          MyImage->FileName, MyImage->M, MyImage->N, NewWidth, NewHeight);

  /* Copy Image Data */
  for(m=0;m<MyImage->M;m++) {
    memcpy(&TempImage->Signal[m+NewStartM][NewStartN*MyImage->ArrayType],
           MyImage->Signal[m],MyImage->N*MyImage->ArrayType*sizeof(float));
  } 

  /* Exchange Signals */
  TempSignal=TempImage->Signal;
  TempImage->Signal=MyImage->Signal;
  MyImage->Signal=TempSignal;

  /* Exchange coords. */
  Swap(MyImage->M,TempImage->M);
  Swap(MyImage->N,TempImage->N);
  
  FreeImage(TempImage);
}


/***************************************************************************
[NAME]
ShrinkImage

[SYNOPSIS]
void ShrinkImage(Image *MyImage,
                 int NewWitdth,
                 int NewHeigth,
                 int NewArea);

[DESCRIPTION]

Shrinks the image to new dimensions {\tt NewWidth $\times$
NewHeight}. The switch {\tt NewArea} determines where the sub-image
should be taken from in the old image. {\tt NewArea} can have the
following values

\begin{verbatim}
NewArea = \_UpperLeft,  \_UpperMiddle,  \_UpperRight
          \_MiddleLeft, \_MiddleMiddle, \_MiddleRight
          \_LowerLeft,  \_LowerMiddle,  \_LowerRight
\end{verbatim}
[USAGE]
{\tt ShrinkImage(Test, 98,98,\_MiddleMiddle);}

Shrinks {\tt Test} to dimensions 98 $\times$ 98. The image is cropped
equally in all directions, i.e. the centrum remains in the center. 

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void ShrinkImage(Image *MyImage,int NewWidth ,int NewHeight,int NewArea)
{
  int NewStartM,NewStartN,m;
  Image *TempImage;
  int Temp;
  float **TempSignal;
  

  if(NewWidth>MyImage->M||NewHeight>MyImage->N)
    Error("ShrinkImage: Dimensions must be smaller than before");

  TempImage=NewFloatImage(MyImage->FileName,
                          NewWidth,NewHeight,MyImage->ArrayType);

  switch(NewArea) {
  case _UpperLeft: case _MiddleLeft: case _LowerLeft:
  	NewStartM=0;	
      	break;
  case _UpperRight: case _MiddleRight: case _LowerRight:
	NewStartM=MyImage->M-NewWidth;
	break;
  case _UpperMiddle: case _MiddleMiddle: case _LowerMiddle:
	NewStartM=(MyImage->M-NewWidth)/2;
	if((!(MyImage->M%2))&&(NewWidth%2)) NewStartM++;
	break;
  default: 
	NewStartM=0;
	Error("ShrinkImage: Unknown Area");
  }

  switch(NewArea) {
  case _UpperLeft: case _UpperMiddle: case _UpperRight:
  	NewStartN=MyImage->N-NewHeight;	
      	break;
  case _MiddleLeft: case _MiddleMiddle: case _MiddleRight:
	NewStartN=(MyImage->N-NewHeight)/2;
	if((!(MyImage->N%2))&&(NewHeight%2)) NewStartN++;
	break;
  case _LowerLeft: case _LowerMiddle: case _LowerRight:
	NewStartN=0;
	break;
  default:
        NewStartN=0;
  }
  /* update image real coords. */
  MyImage->Xmin=MyImage->Xmin+NewStartM*MyImage->DeltaX;
  MyImage->Ymin=MyImage->Ymin+NewStartN*MyImage->DeltaY;

  Print(_DDebug,"ShrinkImage: Shrinking `%s' (%dx%d) to new dimensions (%dx%d)\n",
          MyImage->FileName, MyImage->M, MyImage->N, NewWidth, NewHeight);

  /* Copy Image Data */
  for(m=0;m<NewWidth;m++) {
    memcpy(TempImage->Signal[m],
           &MyImage->Signal[NewStartM+m][NewStartN*MyImage->ArrayType],
                            NewHeight*MyImage->ArrayType*sizeof(float));
  } 

  /* Exchange Signals */
  TempSignal=TempImage->Signal;
  TempImage->Signal=MyImage->Signal;
  MyImage->Signal=TempSignal;

  /* Exchange coords. */
  Swap(MyImage->M,TempImage->M);
  Swap(MyImage->N,TempImage->N);
  
  FreeImage(TempImage);
}


/***************************************************************************
[NAME]
ReadDAT

[SYNOPSIS]
Image *ReadDAT(char *FileName);

[DESCRIPTION]

The function reads a picture in Peter Tofts {\tt .dat} (raw float)
format from the file {\tt 'FileName.dat'}. All memory allocation is done
internally. A pointer to the new image is returned. Included for
compatibility reasons. May not appear in future versions.

[USAGE]
{\tt Test = ReadDAT("Test");}

Reads the image {\tt Test.dat} and return it in {\tt Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
Image * ReadDAT(char *FileName)
{
  int type,m,M,N;
  FILE *InFile;
  char InFileName[100];
  Image *NewImage;

  strcpy(InFileName,FileName);
  strcat(InFileName,".dat");

  Print(_DDebug,"ReadDAT: Reading `%s' ",InFileName);
  
  if(!(InFile=fopen(InFileName,"rb")))	
    Error("Error opening file: `%s'",InFileName);
  
  fread(&type,sizeof(int),1,InFile);
  fread(&M,sizeof(int),1,InFile);
  fread(&N,sizeof(int),1,InFile);

  Print(_DNormal,"(%dx%d)\n", M,N);

  NewImage=NewFloatImage(FileName,M,N,_RealArray);

  for(m=0;m<M;m++)
    if((fread(NewImage->Signal[m],sizeof(float),N,InFile))!=N)
      Error("Error reading file: `%s'",InFileName);

  fclose(InFile);
  InitImage(NewImage);
  return NewImage;
}


/***************************************************************************
[NAME]
WriteDAT

[SYNOPSIS]
void WriteDAT(Image *MyImage);

[DESCRIPTION]

The function writes a picture in {\tt .dat} format. The name of the
saved image is dertermined by {\tt MyImage->FileName}. Only included
for compatibility reasons, may not be included in future versions.

[USAGE]
{\tt WriteDAT(Test);}

Write the image {\tt Test} to a {\tt .dat} file .

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void WriteDAT(Image * MyImage)
{
  int m;
  FILE *OutFile;
  char OutFileName[100];

  if (MyImage->ArrayType==_ComplexArray) 
    Error("Do not use DAT format for Complex images, use FIF instead");

  strcpy(OutFileName,MyImage->FileName);
  strcat(OutFileName,".dat");
  if(!(OutFile=fopen(OutFileName,"wb")))
    Error("Error opening file: `%s'",OutFileName);

  fwrite(&MyImage->Type,sizeof(int),1,OutFile);
  fwrite(&MyImage->M,sizeof(int),1,OutFile);
  fwrite(&MyImage->N,sizeof(int),1,OutFile);

  Print(_DNormal,"WriteDAT: Writing `%s' (%dx%d) \n",
          OutFileName,MyImage->M,MyImage->N);

  for(m=0;m<MyImage->M;m++)
    if((fwrite(MyImage->Signal[m],sizeof(float),
	       MyImage->N,OutFile))!=MyImage->N)
    Error("Error writing to file: `%s'",OutFileName);

  fclose(OutFile);
}

/***************************************************************************
[NAME]
ReadFIFHeader

[SYNOPSIS]
Image *ReadFIFHeader(char *FileName);

[DESCRIPTION]

The function reads the header from a picture in {\tt .fif} format (raw
float with header) from the file {\tt FileName.fif}. A pointer to an
{\tt Image} structure with no image data is returned. 

[USAGE]
{\tt Test = ReadFIF("Test");}

Reads the image {\tt Test.fif} and returns the pointer in {\tt Test}.

[NOTE]

The Image structure returned by this function must not be freed by {\tt
 FreeImage} as it contains no signal information --- use {\tt free} instead.

[REVISION] 
Dec. 94, JJJ and PT
***************************************************************************/
Image *ReadFIFHeader(char *FileName)
{
  int DoConvert;
  FILE *InFile;
  char InFileName[100];
  Image *NewImage;

  DoConvert=FALSE;
    
  strcpy(InFileName,FileName);
  strcat(InFileName,".fif");

  Print(_DNormal,"ReadFIFHeader: Reading `%s' \n",InFileName);
  
  if(!(InFile=fopen(InFileName,"rb")))	
    Error("Error opening file: `%s'",InFileName);
  
  if(!(NewImage=(Image *) malloc(sizeof(Image)))) 
    Error( "\nMemory allocation problems (ReadFIFHeader)\n" );

  fread(NewImage,sizeof(Image),1,InFile);

  if(NewImage->FIFIdType!=_00EI) {
    DoConvert=TRUE;
    convert4((char *)&NewImage->FIFIdType);
    if(NewImage->FIFIdType!=_00EI) Error("Error reading FIF (Converting)\n"); 
  }
  
  /* Converting FIF header.*/
  if(DoConvert) {
    convert4((char *)&NewImage->Type);
    convert4((char *)&NewImage->M);
    convert4((char *)&NewImage->N);
    convert4((char *)&NewImage->ArrayType);
    convert4((char *)&NewImage->Xmin);
    convert4((char *)&NewImage->Ymin);
    convert4((char *)&NewImage->DeltaX);
    convert4((char *)&NewImage->DeltaY);
    convert4((char *)&NewImage->SignalMin);
    convert4((char *)&NewImage->SignalMax);
  }

  fclose(InFile);
  return NewImage;
}


/***************************************************************************
[NAME]
ReadFIF

[SYNOPSIS]
Image *ReadFIF(char *FileName);

[DESCRIPTION]

The function reads a picture in {\tt .fif} format (raw float with
header) from the file {\tt 'FileName.fif'}. All memory allocation is
done internally. A pointer to the new image is returned.

[USAGE]
{\tt Test = ReadFIF("Test");}

Reads the image {\tt Test.fif} and return the pointer in {\tt Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
Image * ReadFIF(char *FileName)
{
  int m,n,DoConvert;
  FILE *InFile;
  char InFileName[100];
  Image *NewImage,*TempImage;

  DoConvert=FALSE;
    
  strcpy(InFileName,FileName);
  strcat(InFileName,".fif");

  Print(_DNormal,"ReadFIF: Reading `%s' \n",InFileName);
  
  if(!(InFile=fopen(InFileName,"rb")))	
    Error("Error opening file: `%s'",InFileName);
  
  if(!(NewImage=(Image *) malloc(sizeof(Image)))) 
    Error( "\nMemory allocation problems (ReadFIF)\n" );

  fread(NewImage,sizeof(Image),1,InFile);

  if(NewImage->FIFIdType!=_00EI) {
    DoConvert=TRUE;
    convert4((char *)&NewImage->FIFIdType);
    if(NewImage->FIFIdType!=_00EI) Error("Error reading FIF (Converting)\n"); 
  }
  
  /* Converting FIF header.*/
  if(DoConvert) {
    convert4((char *)&NewImage->Type);
    convert4((char *)&NewImage->M);
    convert4((char *)&NewImage->N);
    convert4((char *)&NewImage->ArrayType);
    convert4((char *)&NewImage->Xmin);
    convert4((char *)&NewImage->Ymin);
    convert4((char *)&NewImage->DeltaX);
    convert4((char *)&NewImage->DeltaY);
    convert4((char *)&NewImage->SignalMin);
    convert4((char *)&NewImage->SignalMax);
  }
  /*                       */

  TempImage=NewFloatImage(FileName,NewImage->M,NewImage->N,NewImage->ArrayType);
  NewImage->Signal=TempImage->Signal;
  Free(TempImage);

  for(m=0;m<NewImage->M;m++)
    if((fread(NewImage->Signal[m],sizeof(float),
        NewImage->N*NewImage->ArrayType,InFile))!=NewImage->N*NewImage->ArrayType)
      Error("Error reading file: `%s'",InFileName);
  
  if(DoConvert) 
    for(m=0;m<NewImage->M;m++)
      for(n=0;n<(NewImage->N*NewImage->ArrayType);n++)
	convert4((char *)&NewImage->Signal[m][n]);
  
  
  fclose(InFile);
  return NewImage;
}


/***************************************************************************
[NAME]
WriteFIF

[SYNOPSIS]
void WriteFIF(Image *MyImage);

[DESCRIPTION]

The function writes a picture in {\tt .fif} format. The name of the
saved image is determined by {\tt MyImage->FileName}. The file
extension is appended to the filenames in all the `write' routines.

[USAGE]
{\tt WriteFIF(Test);}

Write the image {\tt Test} to a {\tt .fif} file .

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void WriteFIF(Image * MyImage)
{
  int m;
  FILE *OutFile;
  char OutFileName[100];

  strcpy(OutFileName,MyImage->FileName);
  strcat(OutFileName,".fif");

  Print(_DNormal,"WriteFIF: Writing `%s' (%dx%d) \n",
          OutFileName,MyImage->M,MyImage->N);

  if(!(OutFile=fopen(OutFileName,"wb"))) 
    Error("Error opening file: `%s'",OutFileName);

  fwrite(MyImage,sizeof(Image),1,OutFile);

  for(m=0;m<MyImage->M;m++) 
    if((fwrite(MyImage->Signal[m],sizeof(float),
        MyImage->N*MyImage->ArrayType,OutFile))!=MyImage->N*MyImage->ArrayType)
      Error("Error writing to file: `%s'",OutFileName);

  fclose(OutFile);
}


/***************************************************************************
[NAME]
ReadGif

[SYNOPSIS]
Image *ReadGif(char *FileName);

[DESCRIPTION]

The function reads a picture in {\tt .gif}, `Graphics Interchange File
format'. All memory allocation is done internally. A pointer is
returned to the new image.

[USAGE]

{\tt Test = ReadGif("Test");}

Read the image {\tt Test.gif} and return it in {\tt Test}.

[NOTE] 

This function makes use the public domain library `GIFLib' to do the
actual decoding of the images

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
Image * ReadGif(char * FileName)
{
  int i,j,M,N;
  char InFileName[100];
  Image *NewImage;
  GifRecordType Dummy;
  GifPixelType *GifLine;
  GifFileType *GifFile;

  /* Open stdin for the input gif file */
  strcpy(InFileName,FileName);
  strcat(InFileName,".gif");
  if (!(GifFile = DGifOpenFileName(InFileName)))
    Error("Error opening file: `%s'",InFileName);

  Print(_DNormal,"ReadGif: Reading `%s' \n",InFileName);

  /* Is the Image descriptor next ? */
  if (DGifGetRecordType(GifFile,&Dummy) == GIF_ERROR) 
    Error("Error reading file: `%s'",InFileName);

  /* If it is -> read it, else exit */
  if (Dummy==IMAGE_DESC_RECORD_TYPE)
    DGifGetImageDesc(GifFile);
  else 
    Error("Error reading file: `%s'",InFileName);

  M=GifFile->IWidth;
  N=GifFile->IHeight;

  Print(_DNormal,"(%dx%d, %d colors)\n",M,N,(1<<GifFile->IBitsPerPixel));

  /* Allocate memory for one line of data */
  if (!(GifLine=(GifPixelType *)malloc(M*sizeof(GifPixelType))))
    Error("Memory allocation problems (ReadGif)");

  /* Allocate memory for floating point image */
  NewImage = NewFloatImage(FileName,M,N,_RealArray);

  /* Read picture linewise */
  for (i=N-1; i>=0; i--) { 
    if(DGifGetLine(GifFile, &GifLine[0], M) == GIF_ERROR)
      Error("Error reading file: `%s'",InFileName);
    for (j=0; j<M; j++) 
      NewImage->Signal[j][i]=GifLine[j];
  }	

  /* Close infile */
  if (DGifCloseFile(GifFile) == GIF_ERROR)
    Error("Error closing file: `%s'",InFileName);

  Free(GifLine);
  InitImage(NewImage);
  
  return(NewImage);
}


/***************************************************************************
[NAME]
WriteGif

[SYNOPSIS]
void WriteGif(Image *MyImage, 
              char ColorFileName,
              unsigned char ColorDepth,
              float reqmin,
              float reqmax,
              int SaveAsABS);

[DESCRIPTION]

The function writes a picture in the {\tt .gif} format. The number of
colors is specified by {\tt 1<<ColorDepth}.  Picture values are
converted liniarly to match the choosen number of discrete leves in a
{\tt .gif} file.  The palette mapping will be taken from the file {\tt
ColorFileName.map} if it exists, else a generic linear grey scale
palette will be used. The name of the saved image is dertermined by
{\tt MyImage->FileName}. The \tc{gif} format can only be used for real
images. If the image is is complex, the saved image will be a real
(absolute) copy of it. If the two values \tc{reqmin} and \tc{reqmax}
are equal, the picture output values will be between the largest and
smallest value in the picture, else the output values are truncated
into the range \tc{[reqmin;reqmax]}. \tc{SaveAsABS} determines wether
the image should be converted to absolute values. \tc{SaveAsABS}=0
means no converting, \tc{SaveAsABS}=1 converts.

The palette file is an ordinary ASCII file, that contains a line describing
the color for each discrete level in the file, i.e.  if
\tc{ColorDepth} equals 8, there must be 256 lines in the palette file,
each containing three values to describe the color in terms of the
three primary colors, red, green and blue. the values should be
between 0 and 255.  A line to describe the color green should contain
\tc{0 255 0}

[USAGE]
{\tt WriteGif(Test, "gamma1.map", 8);}

Write the image {\tt Test} to a {\tt .gif} file with 256 colors and
using the colormap {\tt gamma1.map}.

[NOTE] 

This function makes use the public domain library `GIFLib' to do the
actual encodin of the images

[REVISION]
Oct. 94, JJJ
***************************************************************************/
void WriteGif(Image *MyImage, char *ColorFileName, unsigned char ColorDepth,
	      float reqmin, float reqmax, int SaveAsABS)
{
  int i,j, Red, Green, Blue, RestoreImage;
  float MaxVal,MinVal,Diff,tempvalue;
  char OutFileName[100];
  short int ColorLevels;
  GifPixelType *GifLine;
  GifColorType *ColorMap;
  GifFileType *GifFile;
  Image *TempImage;
  FILE *ColorFile = NULL;

  RestoreImage=FALSE;

  if ((SaveAsABS==TRUE) || (MyImage->ArrayType==_ComplexArray)) {
    TempImage=CopyImage(MyImage);
    AbsoluteImage(TempImage);
    RestoreImage=TRUE;
  }
  else TempImage=MyImage;

  if (!(GifLine=(GifPixelType *)malloc(TempImage->M*sizeof(GifPixelType))))
    Error("Memory allocation problems (WriteGif)");

  /* Allocate momory for palette */
  ColorLevels=(1 << ColorDepth);
    if ((ColorMap = (GifColorType *) 
                     malloc(sizeof(GifColorType)*ColorLevels)) == NULL)
      Error("Memory allocation problems (WriteGif)");

  /* Read Palette from file */
  if ((ColorFile=fopen(ColorFileName, "rt"))) {
    for (i = 0; i<(1<<ColorDepth); i++) {
      if (feof(ColorFile))
        Error("Not enough colors in palettefile ");
      if (fscanf(ColorFile, "%d %d %d\n", &Red, &Green, &Blue)!=3)
        Error("Error reading file: `%s'",ColorFileName);
      ColorMap[i].Red = Red;
      ColorMap[i].Green = Green;
      ColorMap[i].Blue = Blue;
    }
  }

  /* No palettefile -> initialize greyscale palette */
  else {
  if (ColorFileName[0]!='\0')
    Print(_DNormal, "No color file found, using default grey palette\n");
  for (i = 0, j = 0 ; i < (1<<ColorDepth) ; i++, j+=(256>>ColorDepth)) 
    ColorMap[i].Red=ColorMap[i].Green=ColorMap[i].Blue=j;
  }

  /* Open stdout for the output file: */
  strcpy(OutFileName,TempImage->FileName);
  strcat(OutFileName,".gif");
  if (!(GifFile = EGifOpenFileName(OutFileName, FALSE)))
    Error("Error opening file: `%s'",OutFileName);

  /* Dump out the screen descriptor: */
  if (EGifPutScreenDesc(GifFile, TempImage->M, 
			TempImage->N, ColorDepth, 0, ColorDepth, ColorMap)== GIF_ERROR)
    Error("Error writing screen descriptor to file: `%s'",OutFileName);
  
  /* Dump out the image descriptor: */
  if (EGifPutImageDesc(GifFile,
		       0, 0, TempImage->M, TempImage->N, FALSE, ColorDepth, FALSE) == GIF_ERROR)
    Error("Error writing image descriptor to file: `%s'",OutFileName);
  
  /* Check colorrange in picture */
  if (reqmin==reqmax){
    MaxVal=(MinVal=TempImage->Signal[0][0]);
    for (i=0; i<TempImage->M; i++) 
      for (j=0; j<TempImage->N; j++) {
	if (MinVal>TempImage->Signal[i][j]) MinVal=TempImage->Signal[i][j]; 
	if (MaxVal<TempImage->Signal[i][j]) MaxVal=TempImage->Signal[i][j]; 
      }
  }
  else {
    MinVal=reqmin;
    MaxVal=reqmax;
  }
  
  Print(_DNormal,"WriteGIF: Writing `%s'\n",OutFileName);
  Print(_DDebug,"WriteGIF: (%dx%d, %d colors), Values: %3.2f->%3.2f \n",
	TempImage->M,TempImage->N,(1<<ColorDepth),MinVal,MaxVal);
  
  /* Dump picture linewise */
  Diff=(MaxVal-MinVal)/((1<<ColorDepth)-1);

  for (i=TempImage->N-1; i>=0; i--) {
    for (j=0; j<TempImage->M; j++){
      tempvalue=min( (max( (TempImage->Signal[j][i]) , (MinVal) ) ), (MaxVal) );
      GifLine[j]=(unsigned char)((tempvalue-MinVal)/Diff);
    }
    if(EGifPutLine(GifFile, &GifLine[0], TempImage->M) == GIF_ERROR)
      Error("Error writing lines to file `%s' ",OutFileName);
  }	

  if (EGifCloseFile(GifFile) == GIF_ERROR)
    Error("Error closing file `%s'",OutFileName);

  if (RestoreImage) 
    Free(TempImage);
  
  Free(GifLine);
  Free(ColorMap);
}


/***************************************************************************
[NAME]
ReadMatLab

[SYNOPSIS]
Image *ReadMatLab(char *FileName);

[DESCRIPTION]

The function reads a picture  in MatLab {\tt .mat} format. All memory
allocation  is  done internally. A  pointer  is returned  to the  new
image. The function reads pictures written from DOS, Linux and HP9000.

[USAGE]

{\tt Test = ReadMatLab("Test");}

Read the image {\tt Test.mat} and return it in {\tt Test}.

[REVISION]
Oct. 94, JJJ and PAP \\
Dec. 94, PAP, Last changes

***************************************************************************/

Image * ReadMatLab(char *FileName)
{
  FILE *InFile;
  Fmatrix Header;
  Image * NewImage;
  int m,n,M,N,DoConvert,NumberOfBytes,ElementType;
  char InFileName[100],MatName[100],TempString[15];
  unsigned char InChars[8];

  strcpy(InFileName,FileName);
  strcat(InFileName,".mat");

  if(!(InFile=fopen(InFileName,"rb"))) 
    Error("Error opening file: `%s'",InFileName);
  if(fread(&Header, sizeof(Fmatrix), 1, InFile)!=1)
    Error("Error reading file: `%s'",InFileName);
 
  Print(_DNormal,"ReadMatlab: Reading `%s' ",InFileName);

  DoConvert=FALSE;
  if(((unsigned)Header.type>0x0000ffff)||((Header.type==0)&&(MATLABTYPE==1000))) {
    convert4((char *)&Header.type);
    convert4((char *)&Header.mrows);
    convert4((char *)&Header.ncols);
    convert4((char *)&Header.imagf);
    convert4((char *)&Header.namlen);
    DoConvert=TRUE;
  }
  
  sprintf(TempString,"%04i",Header.type);
  if(strlen(TempString)>4) Error("Wrong Type in Matlab file:%s\n",TempString);
  
  switch(TempString[0]) {
    case '0': 
    if(MATLABTYPE==1000 && DoConvert==FALSE) Error("Wrong type of elements");;
    break;
    case '1':
    if(MATLABTYPE==0 && DoConvert==FALSE) Error("Wrong type of elements");;
    break;
    default:
    Error("This is NOT a VAX or CRAY\n");
    break;
  }
   
  if(TempString[1]!='0') Error("Wrong Type in Matlab file:%s\n",TempString);
  
  ElementType=TempString[2]-'0';
  switch(TempString[2]) {
  case '0': 
    NumberOfBytes=8;
    break;
  case '1':
    NumberOfBytes=4;
    break;
  case '2':
    NumberOfBytes=4;
    break;
  case '3':
    NumberOfBytes=2;
    break;
  case '4':
    NumberOfBytes=2;
    break;
  case '5':
    NumberOfBytes=1;
    break;
  default:
    NumberOfBytes=1;
    Error("Wrong Type in Matlab file:%s\n",TempString);
    break;
  }

  if(TempString[3]!='0') Error("This Matlab type can not be used in this program.\n");

  if(fread(MatName,sizeof(char),Header.namlen,InFile)!=Header.namlen)
    Error("Error reading file: `%s'",InFileName);

  N=Header.ncols;
  M=Header.mrows;

  if(Header.imagf==1) NewImage=NewFloatImage(FileName,M,N,_ComplexArray);
  else NewImage=NewFloatImage(FileName,M,N,_RealArray);

  for (n=0;n<N;n++)
    for (m=0;m<M;m++) {
      if(fread(&InChars,sizeof(char),NumberOfBytes,InFile)!=NumberOfBytes)
        Error("Error reading file: `%s'",InFileName);
      switch(ElementType) {
      case 0: 
	if(DoConvert) convert8(InChars);
	NewImage->Signal[m][n*NewImage->ArrayType]=* (double *)InChars;
	break;
      case 1:
	if(DoConvert) convert4(InChars);
	NewImage->Signal[m][n*NewImage->ArrayType]=* (float *)InChars;
	break;
      case 2:
	if(DoConvert) convert4(InChars);
	NewImage->Signal[m][n*NewImage->ArrayType]=* (int*)InChars;
	break;
      case 3:
	if(DoConvert) convert2(InChars);
	NewImage->Signal[m][n*NewImage->ArrayType]=* (short int *)InChars;
	break;
      case 4:
	if(DoConvert) convert2(InChars);
	NewImage->Signal[m][n*NewImage->ArrayType]=* (unsigned short int *)InChars;
	break;
      case 5:
	NewImage->Signal[m][n*NewImage->ArrayType]=* (unsigned char *)InChars;
	break;
      }
     }

  if(NewImage->ArrayType==_ComplexArray)
    for (n=0;n<N;n++)
      for (m=0;m<M;m++) {
	if(fread(&InChars,sizeof(char),NumberOfBytes,InFile)!=NumberOfBytes)
	  Error("Error reading file: `%s'",InFileName);
	switch(ElementType) {
	case 0: 
	  if(DoConvert) convert8(InChars);
	  NewImage->Signal[m][n*NewImage->ArrayType]=*(double*)InChars;
	  break;
	case 1:
	  if(DoConvert) convert4(InChars);
	  NewImage->Signal[m][n*NewImage->ArrayType]=*(float*)InChars;
	  break;
	case 2:
	  if(DoConvert) convert4(InChars);
	  NewImage->Signal[m][n*NewImage->ArrayType]=*(int*)InChars;
	  break;
	case 3:
	  if(DoConvert) convert2(InChars);
	  NewImage->Signal[m][n*NewImage->ArrayType]=*(short int*)InChars;
	  break;
	case 4:
	  if(DoConvert) convert2(InChars);
	  NewImage->Signal[m][n*NewImage->ArrayType]=*(unsigned short int*)InChars;
	  break;
	case 5:
	  NewImage->Signal[m][n*NewImage->ArrayType]=*(unsigned char*)InChars;
	  break;
	}
      }
  
  fclose(InFile);
  InitImage(NewImage);
  return(NewImage);
}


/***************************************************************************
[NAME]
WriteMatLab

[SYNOPSIS]
void WriteMatLab(Image *MyImage);

[DESCRIPTION]

The function writes a picture in native MatLab {\tt .mat} format. The
name of the saved image is dertermined by {\tt MyImage->FileName}. The
type is determined by the host/server.

[USAGE]
{\tt WriteMatLab(Test);}

Writes the image {\tt Test} to a {\tt .mat} file .

[REVISION]
Oct. 94, JJJ and PAP\\
Dec. 94, PAP, Last changes
***************************************************************************/
void WriteMatLab(Image * MyImage)
{
  FILE *OutFile;
  Fmatrix Header;
  int m,n;
  double TempDouble;
  char OutFileName[100],mnavn[100],* strprt;

  /* Choosing a regular matlab name. */
  strprt=strrchr(MyImage->FileName,'/');
  if(strprt) strcpy(mnavn,strprt+1);
  else strcpy(mnavn,MyImage->FileName);

  Header.type   = MATLABTYPE;
  Header.mrows  = MyImage->M;
  Header.ncols  = MyImage->N;
  Header.imagf  = (MyImage->ArrayType==_ComplexArray);
  Header.namlen = strlen(mnavn)+1;   /* Rember end of string char ! */

  strcpy(OutFileName,MyImage->FileName);
  strcat(OutFileName,".mat");

  Print(_DNormal, "WriteMAT: Writing `%s' (%dx%d) \n",
          OutFileName,MyImage->M,MyImage->N);

  /* Dump header */
  if(!(OutFile=fopen(OutFileName,"wb"))) 
    Error("Error opening file: `%s'",OutFileName);
  if ((fwrite(&Header, sizeof(Fmatrix), 1, OutFile))!=1)
    Error("Error writing file: `%s'",OutFileName);
  if ((fwrite(&mnavn, sizeof(char), 
       (int)Header.namlen, OutFile))!=(int)Header.namlen)
    Error("Error writing file: `%s'",OutFileName);
   
  /* Dump data Real part */
  for (n=0;n<MyImage->N;n++)
    for (m=0;m<MyImage->M;m++) {   
      TempDouble=(double)MyImage->Signal[m][n*MyImage->ArrayType];
      if ((fwrite(&TempDouble, sizeof(double), 1, OutFile))!=1)
        Error("Error writing file: `%s'",OutFileName);
    }
  /* Dump Data Imag. part if any. */
  if(MyImage->ArrayType==_ComplexArray)
    for (n=0;n<MyImage->N;n++)
      for (m=0;m<MyImage->M;m++) {   
	TempDouble=(double)MyImage->Signal[m][n*MyImage->ArrayType+1];
	if ((fwrite(&TempDouble, sizeof(double), 1, OutFile))!=1)
	  Error("Error writing file: `%s'",OutFileName);
      }

  fclose(OutFile);
}


/***************************************************************************
[NAME]
ReadAnalyze

[SYNOPSIS]
Image *ReadAnalyze(char *FileName,
               int LayerNumber);

[DESCRIPTION] 

The function reads a picture in General Electric `Analyze' format from
{\tt FileName.img} with a header file named {\tt FileName.hdr}, the
function reads the image number \tc{layernumber}.  A pointer is
returned to the new image. An analyze file is specified in the
\tc{.ini} file as file type \tc{.analyze}. If \tc{layernumber} is
negative an additional flip of the image is performed. 

[USAGE]

{\tt Test = ReadAnalyze("Test",1);}

Read the analyze file {\tt Test.img} with header {\tt Test.hdr} and
return it in {\tt Test}. The first sinogram is read.

[REVISION]
Oct. 94, PAP\\
Dec. 9, PAP and PT
***************************************************************************/
Image * ReadAnalyze(char *FileName,int LayerNumber)
{
  int TempInt,TempChar,t,m,n,M,N,ABPP,ImageFlag;
  FILE *InFile;
  char InFileName[100];
  Image *NewImage;
  float ScaleFactor,OffSet;
  dsr AnalyzeHeader;

printf("LA=%d\n",LayerNumber);  
  if (LayerNumber<0) {
    LayerNumber=-LayerNumber;
    ImageFlag=TRUE;
  }
  else
    ImageFlag=FALSE;
  /* Reads header file for information about the picture .*/
  strcpy(InFileName,FileName);
  strcat(InFileName,".hdr");
  if(!(InFile=fopen(InFileName,"rb")))	
    Error("Error opening file: '%s'",InFileName);
  
  fread(&AnalyzeHeader,sizeof(dsr),1,InFile);

  fclose(InFile);

  /* Convert from Unix representation of data to PC (Linux).*/

  convert4((char *)&AnalyzeHeader.hk.sizeof_hdr);
  convert4((char *)&AnalyzeHeader.hk.extents);
  convert2((char *)&AnalyzeHeader.hk.session_error);

  for (m=0;m<7;m++) convert2((char *)&AnalyzeHeader.dime.dim[m]);
  convert2((char *)&AnalyzeHeader.dime.bitpix);
  convert2((char *)&AnalyzeHeader.dime.datatype);
  for (m=0;m<7;m++) convert4((char *)&AnalyzeHeader.dime.pixdim[m]);
  convert4((char *)&AnalyzeHeader.dime.funused8);
  convert4((char *)&AnalyzeHeader.dime.funused9);
  convert4((char *)&AnalyzeHeader.dime.glmax);
  convert4((char *)&AnalyzeHeader.dime.glmin);
  /* End Converting */

  /* Byttet om paa x og y i forhold til analyze. */
  M=AnalyzeHeader.dime.dim[2]; 
  N=AnalyzeHeader.dime.dim[1];
  ScaleFactor=AnalyzeHeader.dime.funused9;
  OffSet=AnalyzeHeader.dime.funused8;
  ABPP=AnalyzeHeader.dime.bitpix/8; /* To find number of bytes per pixel.*/

  strcpy(InFileName,FileName);
  strcat(InFileName,".img");
  if(!(InFile=fopen(InFileName,"rb")))	
    Error("Error opening file: '%s'",InFileName);
  
  NewImage=NewFloatImage(FileName,M,N,_RealArray);
  NewImage->DeltaY=AnalyzeHeader.dime.pixdim[1];
  if(ImageFlag==FALSE)
  {
    NewImage->DeltaX=AnalyzeHeader.dime.pixdim[2]*PI/180;
    NewImage->Xmin=0.0;
    NewImage->Ymin=-(N-1)/2*NewImage->DeltaY;
   }
  else
  {
    NewImage->DeltaX=AnalyzeHeader.dime.pixdim[2];
    NewImage->Ymin=-(N-1)/2*NewImage->DeltaY;
    NewImage->Xmin=-(M-1)/2*NewImage->DeltaX;
  }

  /* NOTE !! The offset is NOT implemented */
  NewImage->SignalMax=AnalyzeHeader.dime.glmax*ScaleFactor;
  NewImage->SignalMin=AnalyzeHeader.dime.glmin*ScaleFactor;

  /* Check Layer number, if to large the latest are used.*/
  if(LayerNumber>AnalyzeHeader.dime.dim[3]) {
    Print(_DNormal,"Error reading analyze file: Layer number to large!!\n");
    LayerNumber=AnalyzeHeader.dime.dim[3];
  }

  /* Find the right Layer */
  if(LayerNumber>1) {
    TempInt=M*N*ABPP*(LayerNumber-1);
    if(fseek(InFile,TempInt,SEEK_SET))
      Error("Error reading file: '%s'",InFileName);
  }
  else
    if (LayerNumber==0)
      Error("Layer number equals zero (ReadAnalyze)");

  if (ImageFlag==FALSE)
  {
    for(m=M-1;m>=0;m--)
      for(n=0;n<N;n++) {
	TempInt=0;
	for(t=ABPP;t>0;t--) {
	  TempInt*=256;
	  if((TempChar=fgetc(InFile))==EOF)
	    Error("Error reading file: '%s'",InFileName);
	  TempInt+=(char unsigned)TempChar;
	}
	/* Right sign of the data. */
	if(AnalyzeHeader.dime.glmin<0)
	  if(ABPP==1) TempInt=(signed char)TempInt;
	  else if(ABPP==2) TempInt=(signed short int)TempInt;
	/*   */
	NewImage->Signal[m][n]=(float)TempInt*ScaleFactor;
    }
  }   
  else
  {
    for(m=0;m<M;m++)
      for(n=0;n<N;n++) {
	TempInt=0;
	for(t=ABPP;t>0;t--) {
	  TempInt*=256;
	  if((TempChar=fgetc(InFile))==EOF)
	    Error("Error reading file: '%s'",InFileName);
	  TempInt+=(char unsigned)TempChar;
	}
	/* Right sign of the data. */
	if(AnalyzeHeader.dime.glmin<0)
	  if(ABPP==1) TempInt=(signed char)TempInt;
	  else if(ABPP==2) TempInt=(signed short int)TempInt;
	/*   */
	NewImage->Signal[m][n]=(float)TempInt*ScaleFactor;
    }
    MirrorImage(NewImage);
  }   
    
  fclose(InFile);
  return NewImage;
}


/***************************************************************************
[NAME]
LogImage

[SYN]
void LogImage(Image *MyImage,
              float T0);

[DESCRIPTION]

Calculates  the function  $t=\log\left(\frac{T_0}{T}\right)$ for every
pixel   in the image.   Used   for  converting transmission  scans  to
attenuation factors. if {\tt T0}  is positive, this  value is used for
$T_0$,  else  the $T_0$  is calculated as   the mean of  the {\tt -T0}
lowest lines in the image. 

[USAGE]
{\tt LogImage(Test,-5);}

Converts the image {\tt Test}, with a $T_0$ value calculated as the
average value of the last 5 lines in the image.

[REVISION]
Oct. 94, PAP
***************************************************************************/
void LogImage(Image * MyImage,float T0)
{
  int m,n,Antal;
  float Sum;
  double LogT0;
  
  if(!T0) Error("You must specify T0: Can't be Zero\n");

  if(T0<0) {
    Sum=0;
    Antal=0;
    for(n=0;n<(int)(-T0);n++)
      for(m=0;m<MyImage->M;m++) {
      Sum+=MyImage->Signal[m][n];
      Antal++;
      }
  T0=Sum/Antal;
  }
  Print(_DDebug,"LogImage: Calculatin log(T0/t) with T0:%f\n",T0);
  LogT0=log(T0);
  for(m=0;m<MyImage->M;m++)
    for(n=0;n<MyImage->N;n++) {
    MyImage->Signal[m][n]=LogT0-log(MyImage->Signal[m][n]);
    }
}


/***************************************************************************
[NAME]
ReadImage

[SYNOPSIS]
Image *ReadImage(char *FileName, 
                 int type);

[DESCRIPTION] 

A generic function to read images in a specified format. The avaliable
formats, specified with {\tt type} are {\tt \_GIF}, {\tt \_DAT}, {\tt
\_MAT}, {\tt \_FIF} and {\tt \_Analyze}. The function uses the
previously described routines to read the pictures.

[USAGE] 
{\tt Test = ReadImage("Test",\_GIF);}

Reads the file {\tt Test.gif} in to {\tt Test}.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/

Image * ReadImage(char *filename,int type)
{
  Image *NewImage;

  switch (type) {
    case _DAT:
      NewImage=ReadDAT(filename);
      break;
    case _FIF:
      NewImage=ReadFIF(filename);
      break;
    case _GIF:
      NewImage=ReadGif(filename);
      break;
     case _MAT:
      NewImage=ReadMatLab(filename);
      break;
     case _Analyze:
      NewImage=ReadAnalyze(filename,IniFile.SliceNumber);
      break;
    default:
      NewImage=NULL;
      Error("Image type not supported");
  }
  return NewImage;
}


/***************************************************************************
[NAME]
WriteImage

[SYNOPSIS]
void WriteImage(Image *MyImage, 
                int type);

[DESCRIPTION] 

A generic function to write images in a specified format. The formats
avaliable are {\tt \_GIF}, {\tt \_DAT}, {\tt \_MAT} and {\tt
\_FIF}. The function uses the previously described routines to write
the pictures.


[USAGE] 
{\tt WriteImage(Test,\_GIF);}

Writes {\tt Test} as a {\tt .gif} file.

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/

void WriteImage(Image * MyImage,int type)
{
  switch (type) {
    case _DAT:
      WriteDAT(MyImage);
      break;
    case _FIF:
      WriteFIF(MyImage);
      break;
    case _GIF:
      WriteGif(MyImage,"",8,0,0,FALSE);
      break;
     case _MAT:
      WriteMatLab(MyImage);
      break;
   default:
      Error("Image type not supported");
  }
}


/***************************************************************************
[NAME]
WriteTrace

[SYNOPSIS]
void WriteTrace(Image *MyImage,
                char *TraceFileName,
                int LineNr,
                int Direction);

[DESCRIPTION]

This function writes a trace taken from  the image, from row or column
{\tt LineNr}. The \tc{direction} is  specified with {\tt \_Horizontal}
or {\tt \_Vertical}. The file is a ASCII-file readable by e.g. \tc{gnuplot}.

[USAGE]
{\tt WriteTrace(Test, "testtrace", 50, \_Vertical);}

Writes a trace taken from column 50 in the image {\tt Test}, to the
file {\tt testtrace.trace}.

[REVISION]
Oct. 94, JJJ
***************************************************************************/

void WriteTrace(Image *MyImage, char* TraceFileName, int LineNr, int Direction)
{
  int m;
  float tempfloat;
  FILE *OutFile;
  char OutFileName[100];

  Print(_DNormal,"WriteTrace: Writing trace '%s.trace' from `%s', %s line %d \n", 
	TraceFileName,MyImage->FileName,
	(Direction==_Horizontal) ? "Horisontal" : "Vertical",LineNr);

  strcpy(OutFileName,TraceFileName);
  strcat(OutFileName,".trace");
  if(!(OutFile=fopen(OutFileName,"wb")))
    Error("Error opening file: `%s'",OutFileName);
  
  fprintf(OutFile,"# %s trace from `%s', %s %d\n",
	  (Direction==_Horizontal) ? "Horisontal" : "Vertical", MyImage->FileName, 
	  (Direction==_Horizontal) ? "row" : "collum", LineNr);
  
  if (Direction==_Horizontal)
    for(m=0;m<MyImage->M;m++) {
      if (MyImage->ArrayType==_RealArray) 
	tempfloat=MyImage->Signal[m][LineNr]; 
      else
	tempfloat=sqrt(sq(MyImage->Signal[m][LineNr*2])
		       +sq(MyImage->Signal[m][LineNr*2+1]));
      fprintf(OutFile,"%f %f\n",MyImage->Xmin+MyImage->DeltaX*(float)m,tempfloat);
    }
  else if (Direction==_Vertical) 
    for(m=0;m<MyImage->N;m++) {
      if (MyImage->ArrayType==_RealArray) 
	tempfloat=MyImage->Signal[LineNr][m]; 
      else {
	tempfloat=sqrt(sq(MyImage->Signal[LineNr][m*2])+
		       sq(MyImage->Signal[LineNr][m*2+1]));
      }
      fprintf(OutFile,"%f %f\n",MyImage->Ymin+MyImage->DeltaY*(float)m,tempfloat);
    }
  fclose(OutFile);
}
