/*******************************
*   Prototypes for imgtools.c  *
*******************************/

typedef struct {
  int FIFIdType;       /* Id used to restore FIF:17737:'\0''\0''E''I'   */
  char FileName[100];    /* Name used for saving/restoring this image	 */
  char Description[80];    /* eg. short description of image (optional)	 */
  char Date[6];        /* Date (DDMMYY) (optional)			 */
  int Type;        /* Picture type (eg. format on disc, domain)	 */
  int M;        /* Width of array				 */
  int N;        /* Height of array in no. numbers (Real/Complex) */
  int ArrayType;    /* Defines number fromat: Complex(2) or Real(1)	 */
  float Xmin;        /* Leftmost coor. in original image   		 */
  float Ymin;        /* Lowest coor. in original image	   	 */
  float DeltaX;        /* Quantisation steps in original image	(X)	 */
  float DeltaY;        /* Quantisation steps in original image	(Y)	 */
  float SignalMin;    /* Lowest signalvalue in array			 */
  float SignalMax;    /* Highest signalvalue in array			 */
  float **Signal;    /* The raw data array	 			 */
} Image;

typedef struct {
  int type;        /* type 			*/
  int mrows;        /* row dimension 		*/
  int ncols;        /* column dimension 		*/
  int imagf;        /* flag indicating imag part 	*/
  int namlen;        /* name length (including NULL)	*/
} Fmatrix;        /* Matlab 			*/

#ifdef __hp9000s700
#define MATLABTYPE 1000
#else
#define MATLABTYPE 0
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef PI
#define PI M_PI
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef PI2
#define PI2 M_PI_2
#endif

#define _00EI 17737

/***** Generic **************/
#define TRUE        1
#define FALSE        0
/****************************/

/***** WriteTrace ***********/
#define _Horizontal     1
#define _Vertical       0
/****************************/

/**** Image Tools ***********/
#define _MAT        101     /* Types used in image format. */
#define _GIF        103
#define _DAT        104
#define _FIF        105
#define _Analyze        106
#define _RealArray    1
#define _ComplexArray   2
/****************************/


/******** GetArg ************/
#define _ArgError    0
#define _NoError    1
/****************************/


/*** Shrink/StretchImage ****/
#define _UpperLeft    1
#define _UpperMiddle    2
#define _UpperRight    3
#define _MiddleLeft    4
#define _MiddleMiddle    5
#define _MiddleRight    6
#define _LowerLeft    7
#define _LowerMiddle    8
#define _LowerRight    9
/****************************/

#define Swap(a, b) {Temp=(a);(a)=(b);(b)=Temp;}
#define Free(x) {if(x) free(x); else Error("Trying to free NULL pointer");}

extern void ImgError(char *);

extern Image *NewFloatImage(char *, int, int, int);
extern void FreeImage(Image *);
extern void ZeroImage(Image *);
extern void NormImage(Image *, float, float);
extern void Real2ComplexImage(Image *MyImage);
extern void AbsoluteImage(Image *);
extern void RealImage(Image *);
extern void ImagImage(Image *);
extern void RenameImage(Image *, char *);
extern void MirrorImage(Image *);
extern void StretchImage(Image *, int, int, int);
extern void ShrinkImage(Image *, int, int, int);
extern Image *CopyImage(Image *);
extern void CropImage(Image *, int, int, int, int);

extern Image *ReadDAT(char *);
extern Image *ReadFIFHeader(char *);
extern Image *ReadFIF(char *);
extern Image *ReadGif(char *);
extern Image *ReadMatLab(char *);
extern Image *ReadAanalyze(char *FileName, int LayerNumber);
extern Image *ReadImage(char *, int);

extern void WriteDAT(Image *);
extern void WriteFIF(Image *);
extern void WriteGif(Image *, char *, unsigned char, float, float, int);
extern void WriteMatLab(Image *);
extern void WriteImage(Image *, int);
extern void WriteTrace(Image *, char *, int, int);

extern void LogImage(Image *, float);










