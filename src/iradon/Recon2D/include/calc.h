/***************************
*   Prototypes for calc.c  *
***************************/

#define myPI 3.14159265358979

#define _FFT           -1
#define _IFFT        1

#define _Ramp           100
#define _Hanning        110
#define _Hamming        120

#define _Median         10
#define _Middle         11

typedef float NumberType;

extern void Forward(Image *MyImage, Image *Sinogram);
extern void BackProject(Image *Sinogram, Image *InvMyImage);
extern Image *IChirpSpectrum(Image *MyImage);
extern Image *IFFTSpectrum(Image *MyImage);
extern Image *CentralSliceNN(Image *MyImage);
extern Image *CentralSliceBL(Image *MyImage);
extern Image *CentralSliceCZ(Image *MyImage);
extern Image *FilteredBack(Image *MyImage);
extern Image *BackFilter(Image *MyImage);
extern void FilterRealSpatial(Image *MyImage, int Mode);
extern void FilterRealSpectrum(Image *MyImage, int Mode);
extern float *ComplexChirpZ(float *inarr, int N, int M,
                            float theta0, float phi0, int isign);
extern void RealFFT(float *, int, int);
extern void ComplexFFT(float *, int, int);
extern void VerticalFFT(Image *MyImage, int Isign);
extern void RealVerticalFFT(Image *MyImage, int Isign);
extern void ComplexVerticalFFT(Image *MyImage, int Isign);
extern void FFTImage(Image *MyImage, int Isign);
extern void CheckFFTLength(int);
extern void FFTShift(Image *, int);
extern void sort(int, float *);
extern void ImageFiltering(Image *, int, int, int, int);
void ComplexNdimFFT(float *, unsigned long *, int, int);























