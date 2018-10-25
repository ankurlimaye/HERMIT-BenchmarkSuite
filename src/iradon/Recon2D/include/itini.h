/* defines for it.c */

/****** RadonKernels *******/
#define   _SINC      0
#define   _NN        1
#define   _RNN       2
#define   _RL        3
#define   _P1        4
#define   _P2        5
/***************************/

/****** Algorithms *********/
#define   _ART       0
#define   _EM        1
#define   _CG        2
/***************************/

#define MaxSaves 300

typedef struct{
  int XSamples;             /* No of samples on 1. axis of recon image  */
  int YSamples;             /* No of samples on 2. axis of recon image  */
  float DeltaX;             /* Sampling distance 1. axis of recon image */
  float DeltaY;             /* Sampling distance 2. axis of recon image */
  float Xmin;               /* 1. sample position 1.axis of recon image */
  float Ymin;               /* 1. sample position 2.axis of recon image */
  int ThetaSamples;         /* No of angular samples in sinogram        */
  int RhoSamples;           /* No of spatial samples in sinogram        */
  float DeltaRho;           /* Spatial sampling distance in sinogram    */
  float DeltaTheta;         /* Angular sampling distance in sinogram    */
  float ThetaMin;           /* 1. angular position in sinogram          */
  float RhoMin;             /* 1. spatial position in sinogram          */
  float LowestALevel;       /* Lowest level allowed in A Matrix (0<x<1) */
  int   RadonKernel;        /* How to calculate the A Matrix            */
                            /* 'SINC', 'NN', 'RNN', 'RL', 'P1', 'P2'    */
  int OverSamp;      /* Use oversampling - Using squared no. of samples */
  float Regularization;     /* Use regularisation if not zero.          */
                            /* If positive use Laplace regularization   */
                            /* If negative use identity matrix          */
  int Iterations;           /* Number of Iterations                     */
  int Algorithm;            /* Inversion algorithm : 'EM', 'ART'        */
  int IsFast;               /* Should the fsst/slow algorithm be used   */
  float ConstrainMin;       /* Lower reconstructed constrain limit      */
  float ConstrainMax;       /* Upper reconstructed constrain limit      */
  float Alpha;              /* Algorithm dependent parameter 1          */
  float Beta;               /* Algorithm dependent parameter 2          */
  int IterationType;        /* Algorithm dependent parameter            */
  int SaveMatLab; /* Save the A-Matrix as a MATLab sparse matrix (0|1)  */
  char StartFileName[200];  /* Startguess (optional)                    */
  char RefFileName[200];    /* Analytical constructed image (optional)  */
  char KernelFileName[200]; /* A-Matrix name                            */
  char InFileName[200];     /* Sinogram to be reconstructed             */
  char OutFileName[200];    /* Filename for reconstructed image         */
  int SaveIterations;       /* Number of iterations pr. save 1-200      */
} itINItype;










