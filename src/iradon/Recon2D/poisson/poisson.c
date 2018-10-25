/*********************************************************************
[HEADER]
The program {\sf poisson.c} will simulate a PET scanning.
The program will require a (analytical) sinogram as input.
The program will generate a set of sinograms as function of time.
The simulated sinograms will be Poisson distributed as in real life.
The sinogram is normalized with the value {\tt Norm}.

21/12-94 Peter Toft
*********************************************************************/
#include "iradon.h"

Image *Sim_Sinogram,*Org_Sinogram;
char  Nameorg[80],Name[80],Buffer[50];
float **poi,**Lambda;
int *idum;

extern float poidev(float,int*); /* Numerical Recipes function */
extern float ran1(int*);         /* Numerical Recipes function */
extern float gammln(float);      /* Numerical Recipes function */


/********************************************************************
[NAME]
Generate

[SYNOPSIS]
void Generate(char *SinogramFilename,float Norm,int I)

[DESCRIPTION]
The function {\tt Generate} will generate a set of {\tt I} images,
{\tt Sim1},... In every pixel the increment of the pixel is corrupted 
by Poisson distributed noise. 

[REVISION]
21/12-94 Peter Toft 
********************************************************************/
void Generate(char *SinogramFileName,float Norm)
{
  int r=0,t,R,T;
  float lambda,sum,average,rel;
  char *strptr;
  float *p1,*p2;

  if ((strptr=strrchr(SinogramFileName,'.'))!=NULL)
    SinogramFileName[strptr-SinogramFileName]='\0';
 
  Org_Sinogram=ReadFIF(SinogramFileName);
  T=Org_Sinogram->M;
  R=Org_Sinogram->N;
  Lambda=Org_Sinogram->Signal;
  
  sprintf(Nameorg,"%s.%.0f",SinogramFileName,Norm);
  Sim_Sinogram=CopyImage(Org_Sinogram);
  ZeroImage(Sim_Sinogram);
  poi=Sim_Sinogram->Signal;
  rel=1e2/T;
  printf("Normalizing Image (%i*%i) to the desired average = %f\n",T,R,Norm);
  sum=0;
  for (t=0;t<T;t++)
  {
    p1=Lambda[t];
    for (r=0;r<R;r++)
    {
      lambda=p1[r];
      if (lambda<0) 
        p1[r]=0.0;
      else
        sum+=lambda;
    }
  }

  average=sum/((float)T*R);

  for (t=0;t<T;t++)
  {
    p1=Lambda[t];
    for (r=0;r<R;r++)
      p1[r]*=Norm/average;
  }
  
  RenameImage(Sim_Sinogram,Nameorg);
  for (t=0;t<T;t++) {
    Print(_DNoLog,"%6.2f%% Done\r",t*rel);
    p1=Lambda[t];
    p2=poi[t];
    for (r=0;r<R;r++)
      p2[r]=(int)poidev(p1[r],idum);
  }
  Print(_DNormal,"%6.2f%% Done\n",t*rel);
  WriteFIF(Sim_Sinogram);
}  

void main(int argc,char * argv[])
{
  idum=(int *)malloc(sizeof(int));
  idum[0]=37;
  DebugNiveau=_DNormal;
  strcpy(LogFileName,argv[0]);
  OpenLog(LogFileName);
  if (argc!=3)
  {
    printf("Number of arguments detected %i\n",argc);
    
    puts("Usage -> Poisson Sinogramfile Average");
    exit(1);
  }
  Generate(argv[1],atof(argv[2]));
}


