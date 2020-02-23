/*********************************************************************
[HEADER]

RadonAna version 1.1

Program used to calculate synthetic images and their
corresponding Radon transform. The program uses simple
figures as base elements placed in a $(x,y)$-coordinate
system.

From the author a program {\tt RadonAna} is available which have been
used to generate Radon transform pairs. The basic idea is to specify a
set of scaled, rotated, and shifted primitives from which the Radon
transform is calculated and sampled in both the image domain and the
Radon domain (sinogram).

The program uses the following parameters
\begin{description}
\item[{\tt OutFileName}] [String] With this parameter the base file name is
determined. If {\tt OutFileName} is {\tt man.fif} then the image is
written in {\tt mana.fif} and the sinogram in the file {\tt manr.fif}.
\item[{\tt XSamples}] [Integer] Specifies the number of image samples on the
first axis. Should be an odd number.
\item[{\tt YSamples}] [Integer] Specifies the number of image samples on the
second axis. Should be an odd number. If not given the program uses
{\tt XSamples=YSamples}.
\item[{\tt DeltaX}] [Float] Specifies the sampling distance of both axes in
the image.
\item[{\tt Xmin}] [Float] Specifies the minimum sample position in the image
on the first axis. If not given the image is centered around the middle.
\item[{\tt Ymin}] Specifies the minimum sample position in the image
on the second axis. If not given the image is centered around the middle.
\item[{\tt RhoSamples}] [Integer] Specifies the number of samples in the
distance parameter $\rho$ in the sinogram. Should be an odd number.
\item[{\tt ThetaSamples}] [Integer] Specifies the number of samples in
the angular parameter $\theta$ in the sinogram. The sinogram is
samples linearly from $0$ to (approximately) $\pi$ radians. The
sampling distance in $\theta$ is $\pi/{\tt ThetaSamples}$.
\item[{\tt DeltaRho}] [Float] Specifies the sampling distance in
$\rho$. The program will center the sampling points around $0$.
\item[{\tt OverSamp}] [Integer] This parameter specifies whether an
oversampled image and sinogram should be generated or not {\tt
OverSamp$=0$}. If, e.g., {\tt OverSamp=2} an image is
calculated with 5 times (2*{\tt OverSamp}+1) the resolution in both
parameters. This image is then averaged so the output value make a better
approximation to the average value under the sample. This technique is
also applied the the sinogram. It will reduce aliasing problems.
\item[{\tt GenerateImage}] [Integer] If set to 1 the image is generated.
\item[{\tt GenerateRadon}] [Integer] If set to 1 the sinogram is generated.
\item[{\tt NumberOfShapes}] [Integer] Here the number of primitives are
specified.
\item[{\tt DebugNiveau}] [String] Log files is generated if this parameter is
set to {\tt \_Debug}. No log-file and no screen output can be used
with {\tt DebugNiveau=\_HardCore}.
\end{description}
The ordering of these parameters is arbitrary. After these parameters
a number of lines must follow specifying a scaling, rotation, and
shifting parameters of each primitive. The line uses the following
parameters in this order
\begin{description}
\item[{\tt Shape}$i$] [String] where $i$ is 0 to {\tt
NumberOfShapes}-1.
\item[{\tt Type}] [Integer] The type of primitive, where
\begin{description}
\item[{\tt 1}] A circle centered around (0,0) with radius 1 and
     uniform signal 1 on the disk.
\item[{\tt 2}] A square centered around (0,0). $-1<x<1$ and $-1<y<1$.
     The signal is 0 on the disk.
\item[{\tt 3}] A Gaussian bell $\exp(-x^2-y^2)$.
\item[{\tt 4}] An even sided triangle centered around $(0,0)$. The
     triangle has corners $(1,0)$, $(-\frac{1}{2},\frac{\sqrt{3}}{2})$
     and $(-\frac{1}{2},-\frac{\sqrt{3}}{2})$.
\item[{\tt 5}] The image and the Radon domain is con-terminated with
     white Gaussian noise. Note that in this case the noise term 
     in the Radon domain is NOT the Radon transform of the noise
     in the image.
\item[{\tt 6}] A "pyramid" $f(x,y)=(1-|x|)(1-|y|)$ for $-1<x<1$ and $-1<y<1$.  
\end{description} 
\item[{\tt a}] [Float] Scaling of the primitive in the direction of
the first axis.
\item[{\tt b}] [Float] Scaling of the primitive in the direction of
the second axis.
\item[{\tt x0}] [Float] Shift of the primitive in direction of the first axis.
\item[{\tt y0}] [Float] Shift of the primitive in direction of the second axis.
\item[{\tt phi}] [Float] Rotation of the primitive measured in degrees.
\item[{\tt power}] [Float] Scaling of the intensity of the primitive.
\item[{\tt Description}] [String] Optional description of the primitive.
\end{description}
The program is executed with the initialization file or {\tt .INI}-file 
as sole argument. A typical {\tt .INI}-file could be the
one used to generate the man shown below. 

{\small\begin{verbatim}
OutFileName=man.fif   
XSamples=401          Number of X samples
DeltaX=0.04           Sampling Distance in X
RhoSamples=351        Number of rho samples 
OverSamp=0            Half oversampling factor.   
GenerateImage=1       Generates image if 1. Optional -> Generates if omitted
GenerateRadon=1       Generates sinogram if 1. Optional -> Generates if omitted
DeltaRho=0.05         Sampling distance in rho
ThetaSamples=300      Number of theta samples. Distributed from 0 to pi.
NumberOfShapes=22     Number of simple elements.
DebugNiveau=_DDebug   

        type  a      b     x0    y0      phi power  description
		  
Shape0=  1   1      1.5    0.0   3.7      0   0.8   head
Shape1=  1   0.4    0.2    0.0   2.7      0  -0.2   mouth
Shape2=  2   1.5    2.0    0.0   0.0      0   1.0   stomach
Shape3=  2   1.5    0.3    0.0  -1.0      0  -0.3   belt
Shape4=  2   1.5    0.5    3.3   2.5     30   1.0   right arm  
Shape5=  2   1.5    0.5   -2.4   0.6   -100   1.0   left arm
Shape6=  4   1      1.0    0.0   0.0     90  -0.2   triangle on stomach
Shape7=  1   0.2    0.3   -1.2   3.7      0   0.8   left ear
Shape8=  1   0.2    0.3    1.2   3.7      0   0.8   right ear
Shape9=  2   1.5    0.15   0.0  5.35      0   0.5   lower part of hat
Shape10= 2   0.7    0.6    0.0   6.1      0   0.5   upper part of hat
Shape11= 1   0.3    0.2   -0.5   4.0      0  -0.8   left eye
Shape12= 1   0.3    0.2    0.5   4.0      0  -0.8   right eye
Shape13= 1   0.2    0.4    0.0   3.5      0  -0.4   nose
Shape14= 2   2.0    0.5    2.0  -4.0    -60   1.0   right leg
Shape15= 1   0.5    0.4    5.5   3.7     30   0.5   right hand
Shape16= 1   1      0.5    4.0  -6.5     10   0.5   left foot
Shape17= 2   2.0    0.5   -2.0  -4.0     60   1.0   left leg
Shape18= 1   0.5    0.4   -2.8  -1.6   -100   0.5   left hand
Shape19= 1   1      0.5   -4.0  -6.5    -10   0.5   right foot
Shape20= 4   1      0.5    5.5   4.2    -90   0.5   icecream cone in right hand
Shape21= 1   0.5    0.5    5.5   4.9      0   0.5   only one ice cream cube!
\end{verbatim}}

The program will produce two files {\tt FILENAMEa.fif} and
{\tt FILENAMEr.fif}, which is the image and the radon space.

PT and PAP 25/10-1994\\
PT 14/11-94\\
PT 8/12-94\\
PT 17/2-96\\
PT 12/4-96 More text info above
*********************************************************************/
#include "iradon.h"

extern float gasdev(int *);
extern int *IntVector(int);

struct einfo {
  float a, b, x0, y0, phi, pow;
  int phtype;
};

Image *ImImage, *RtImage, *OutImImage, *OutRtImage;
FILE *indfil;
int XSamples, YSamples, RhoSamples, NumberOfShapes, rmin, ThetaSamples;
int IOkode, OverSamp, GenerateRadon, GenerateImage;
int *idum;
float DeltaX, DeltaRho;
struct einfo *e;
char outname[80];
float sqrt3, topi, topi3, pi3, sqrt3_2, sqrt2;
float Xmin, Ymin;
extern char *IniBuffer;
char OutFileName[200];
char *IniBuffer;
INI IniFile;
char typear[7][50] = {{"Not        "},
                      {"Ellipse    "},
                      {"Square     "},
                      {"Triangle   "},
                      {"Gaus. Bell "},
                      {"Gaus. Noise"},
                      {"Pyramid    "}};

/**********************************************************
[NAME]
image\_func

[SYNOPSIS]
float image_func(float x,float y,struct einfo w)

[DESCRIPTION]
The function will return the image value i the point $(x,y)$
of the chosen base function shifted, rotated and stretched
according to the parameters in {\tt w}.

[NOTE]
Now six functions are supported

[REVISION]
8/12-94 Peter Toft
**********************************************************/
float image_func(int xi, int yi, float x, float y, struct einfo w) {
  float cphi, sphi, xt, yt, base, fx, fy;

  xt = x - w.x0;
  yt = y - w.y0;
  cphi = cos(w.phi);
  sphi = sin(w.phi);
  x = (xt * cphi + yt * sphi) / w.a;
  y = (yt * cphi - xt * sphi) / w.b;

  fx = fabs(x);
  fy = fabs(y);
  switch (w.phtype) {
    case 1 : /* Circular disc */
      if (x * x + y * y < 1.0)
        base = w.pow;
      else
        base = 0.0;
      break;
    case 2 : /* Square */
      if ((fx < 1.0) && (fy < 1.0))
        base = w.pow;
      else
        base = 0.0;
      break;
    case 3 : /* Gaussian bell */
      base = exp(-x * x - y * y) * w.pow;
      break;
    case 4 :
      if ((x < -0.5) || (x > 1.0))
        base = 0.0;
      else if (fy < (1.0 - x) / sqrt3)
        base = w.pow;
      else
        base = 0.0;
      break;
    case 5 : /* White gaussian noise */
      base = gasdev(idum) * w.pow;
      break;
    case 6 : /* "Pyramid" */
      if ((fx < 1.0) && (fy < 1.0))
        base = w.pow * (1.0 - fx) * (1.0 - fy);
      else
        base = 0.0;
      break;
    default :base = 0;
      Error("Function of wrong type detected\n");
  }
  return base;
}

/**********************************************************
[NAME]
radon\_func

[SYNOPSIS]
float radon\_func(float rho,float theta,struct einfo w)

[DESCRIPTION]
The function will return the Radon transform value i the point
$(\rho,\theta)$ of the chosen base function shifted,
rotated and stretched according to the parameters in {\tt w}.

[NOTE]
Now six functions are supported

[REVISION]
8/12-94 Peter Toft 
**********************************************************/
float radon_func(float rho, float theta, struct einfo w) {
  float xp1, xm1, yp1, ct, st, tt, xm, xp, ym, yp;
  float base, gamma, a, b, I1, I2, I3, I4, y1, y0;

  rho -= w.x0 * cos(theta) + w.y0 * sin(theta);
  theta -= w.phi;
  a = w.a;
  b = w.b;
  while (theta >= PI) {
    rho = -rho;
    theta -= PI;
  }
  while (theta < 0) {
    rho = -rho;
    theta += PI;
  }
  ct = cos(theta);
  st = sin(theta);
  gamma = hypot(a * ct, b * st);
  rho = rho / gamma;
  if (theta >= PI2)
    theta = atan(b / a * st / ct) + PI;
  else
    theta = atan(b / a * st / ct);

  switch (w.phtype) {
    case 1 : /* Circular disc */
      if ((tt = 1.0 - rho * rho) > 0.0)
        base = 2.0 * sqrt(tt);
      else
        base = 0.0;
      break;
    case 2 : /* Square */
      if (fabs(rho) > sqrt2)
        base = 0.0;
      else {
        if (rho < 0) {
          rho = -rho;
          theta += PI;
        }
        while (theta >= 2 * PI) theta -= 2 * PI;
        if (theta >= PI) theta = 2 * PI - theta;
        if (theta >= PI * 0.5) theta = PI - theta;
        if (theta >= PI * 0.25) theta = PI * 0.5 - theta; /* Now theta <PI/4 */
        ct = cos(theta);
        st = sin(theta);
        tt = st / ct;
        xp1 = rho / ct - tt;
        xm1 = rho / ct + tt;
        if (theta < 1e-20)
          if (rho < 1)
            base = 2;
          else
            base = 0.0;
        else {
          yp1 = (rho - ct) / st;
          if (xp1 > 1)
            base = 0.0;
          else if (xm1 < 1)
            base = 2 / ct;
          else
            base = hypot(1 - xp1, 1 - yp1);
        }

      }
      break;
    case 3 : /* Gaussian bell */
      base = exp(-rho * rho) * sqrt(PI);
      break;
    case 4 : /* Even side triangle */
      if (rho < 0) {
        rho = -rho;
        theta += M_PI;
      }
      if (rho > 1)
        base = 0.0;
      else {
        while (theta > topi3) theta -= topi3;
        while (theta < 0.0) theta += topi3;
        if (theta > pi3) theta = topi3 - theta;
        ct = cos(theta);
        st = sin(theta);
        if (rho > ct)
          base = 0.0;
        else {
          xp = (rho + st / sqrt3) / (st / sqrt3 + ct);
          yp = (xp - 1.0) / sqrt3;
          if ((rho + 0.5 * ct) > (st * sqrt3_2)) {
            xm = (rho - st / sqrt3) / (ct - st / sqrt3);
            ym = -(xm - 1.0) / sqrt3;
            base = hypot((xp - xm), (yp - ym));
          } else
            base = hypot((xp + 0.5), (yp - (rho + 0.5 * ct) / st));
        }
      }
      break;
    case 5 : /* Noise is added */
      base = gasdev(idum);
      break;
    case 6 : /* "Pyramid" */
      if (fabs(rho) > sqrt2)
        base = 0.0;
      else {
        if (rho < 0)
          rho = -rho;
        while (theta >= 2 * PI) theta -= 2 * PI;
        if (theta >= PI) theta = 2 * PI - theta;
        if (theta >= PI * 0.5) theta = PI - theta;
        if (theta >= PI * 0.25) theta = PI * 0.5 - theta; /* Now theta <PI/4 */
        ct = cos(theta);
        st = sin(theta);
        if (rho > st + ct)
          base = 0.0;
        else if (st < 1e-7)
          base = (1 - rho);
        else {
          if (rho > -st + ct)  /* y1 > -1 */
          {
            y1 = (rho - ct) / st;
            I1 = ct * (1 - y1);
            if (y1 < 0)
              I2 = -0.5 * ct * (1 + y1 * y1);
            else
              I2 = -0.5 * ct * (1 - y1 * y1);
            if (rho > st) {
              I3 = rho * (0.5 * y1 - 1) + 0.5 * (st + y1 * ct);
              if (y1 < 0)
                I4 = 0.5 * rho + (y1 * y1 * (ct + 0.5 * rho) - st) / 3.0;
              else
                I4 = 0.5 * rho - (y1 * y1 * (ct + 0.5 * rho) + st) / 3.0;
            } else {
              y0 = rho / st;
              I3 = rho * (1 - y0 + 0.5 * y1) + 0.5 * (y1 * ct - st);
              if (y1 < 0)
                I4 = (st + y0 * y0 * rho + y1 * y1 * (ct + 0.5 * rho)) / 3 - 0.5 * rho;
              else
                I4 = (st + y0 * y0 * rho - y1 * y1 * (ct + 0.5 * rho)) / 3 - 0.5 * rho;
            }
          } else  /* y1 < -1 */
          {
            I1 = 2 * ct;
            I2 = -ct;
            if (rho > st) {
              I3 = -2 * rho;
              I4 = rho;
            } else {
              y0 = rho / st;
              I3 = -st - rho * rho / st;
              I4 = (2 * st + rho * y0 * y0) / 3;
            }
          }
          base = (I1 + I2 + I3 + I4) / (ct * ct);
        }
      }
      break;
    default :Print(_DNormal, "Function of wrong type detected\n");
      exit(1);
  }

  return base * w.pow * a * b / gamma;
}

/**********************************************************
[NAME]
init

[SYNOPSIS]
void init(void)

[DESCRIPTION]
The function will initialize and read the {\tt .ini}-file.

[REVISION]
25/10-94 Peter Toft 
**********************************************************/
void init(void) {
  int i;
  char Temp[100], Tempi[100];
  char *strptr;

  sqrt3 = sqrt(3.0);
  topi = 2 * PI;
  topi3 = topi / 3.0;
  pi3 = PI / 3.0;
  sqrt3_2 = sqrt3 * 0.5;

  if (!(GetArg(IniBuffer, "DebugLevel", Temp))) DebugNiveau = _DNormal;
  else if (strstr(Temp, "Normal")) DebugNiveau = _DNormal;
  else if (strstr(Temp, "Debug")) DebugNiveau = _DDebug;
  else if (strstr(Temp, "NoScreen")) DebugNiveau = _DNoScreen;
  else if (strstr(Temp, "NoLog")) DebugNiveau = _DNoLog;
  else if (strstr(Temp, "HardCore")) DebugNiveau = _DHardCore;
  else Error("Unknown DebugLevel: '%s'", Temp);

  strcpy(LogFileName, "radonana");

  if (!(GetArg(IniBuffer, "OutFileName", Temp)))
    Error("'OutFileName' entry missing in INI file");
  sscanf(Temp, "%s", OutFileName);

  if (!(strptr = strrchr(OutFileName, '.')))
    Error("Error in IniFile: '%s'", OutFileName);
  OutFileName[strptr - OutFileName] = '\0';
  Print(_DNormal, "'OutFileName' entry is set to %s\n", OutFileName);

  if (!(GetArg(IniBuffer, "XSamples", Temp)))
    Error("'XSamples' entry missing in INI file");
  sscanf(Temp, "%i", &XSamples);
  Print(_DNormal, "'XSamples' entry is set to %i\n", XSamples);

  if (!(GetArg(IniBuffer, "YSamples", Temp))) {
    YSamples = XSamples;
    Print(_DNormal, "'YSamples' entry missing in INI file. Setting it to %i\n", YSamples);
  } else {
    sscanf(Temp, "%i", &YSamples);
    Print(_DNormal, "'YSamples' entry is set to %i\n", YSamples);
  }

  if (!(GetArg(IniBuffer, "DeltaX", Temp)))
    Error("'DeltaX' entry missing in INI file");
  sscanf(Temp, "%f", &DeltaX);
  Print(_DNormal, "'DeltaX' entry is set to %f\n", DeltaX);

  if (!(GetArg(IniBuffer, "Xmin", Temp))) {
    Xmin = (1 - XSamples) * 0.5 * DeltaX;
    Print(_DNormal, "'Xmin' entry missing in INI file. Setting it to %f\n", Xmin);
  } else {
    sscanf(Temp, "%f", &Xmin);
    Print(_DNormal, "'Xmin' entry is set to %f\n", Xmin);
  }

  if (!(GetArg(IniBuffer, "Ymin", Temp))) {
    Ymin = (1 - YSamples) * 0.5 * DeltaX;
    Print(_DNormal, "'Ymin' entry missing in INI file. Setting it to %f\n", Ymin);
  } else {
    sscanf(Temp, "%f", &Ymin);
    Print(_DNormal, "'Ymin' entry is set to %f\n", Ymin);
  }

  if (!(GetArg(IniBuffer, "RhoSamples", Temp)))
    Error("'RhoSamples' entry missing in INI file");
  sscanf(Temp, "%i", &RhoSamples);
  Print(_DNormal, "'RhoSamples' entry is set to %i\n", RhoSamples);

  if (!(GetArg(IniBuffer, "ThetaSamples", Temp)))
    Error("'ThetaSamples' entry missing in INI file");
  sscanf(Temp, "%i", &ThetaSamples);
  Print(_DNormal, "'ThetaSamples' entry is set to %i\n", ThetaSamples);

  if (!(GetArg(IniBuffer, "DeltaRho", Temp)))
    Error("'DeltaRho' entry missing in INI file");
  sscanf(Temp, "%f", &DeltaRho);
  Print(_DNormal, "'DeltaRho' entry is set to %f\n", DeltaRho);

  if (!(GetArg(IniBuffer, "OverSamp", Temp))) {
    Print(_DNormal, "'OverSamp' entry missing in INI file. Setting it to 0\n");
    OverSamp = 0;
  } else {
    sscanf(Temp, "%i", &OverSamp);
    Print(_DNormal, "'OverSamp' entry is set to %i\n", OverSamp);
  }

  if (!(GetArg(IniBuffer, "NumberOfShapes", Temp)))
    Error("'NumberOfShapes' entry missing in INI file");
  sscanf(Temp, "%i", &NumberOfShapes);
  Print(_DNormal, "'NumberOfShapes' entry is set to %i\n", NumberOfShapes);

  e = (struct einfo *) malloc(NumberOfShapes * sizeof(struct einfo));

  for (i = 0 ; i < NumberOfShapes ; i++) {
    sprintf(Tempi, "Shape%i", i);
    if (!(GetArg(IniBuffer, Tempi, Temp))) {
      Print(_DNormal, "Shape Number %i is missing in the INI file. Specify '%s'\n", i, Tempi);
      Error("Terminating Program\n");
    }
    sscanf(Temp, "%i %f %f %f %f %f %f",
           &e[i].phtype, &e[i].a, &e[i].b, &e[i].x0,
           &e[i].y0, &e[i].phi, &e[i].pow);

    e[i].phi = e[i].phi / 180 * PI;
    Print(_DNormal, "%2i : type=%s a=%7.3f b=%7.3f x0=%7.3f y0=%7.3f  pow=%7.3f\n",
          i, typear[e[i].phtype], e[i].a, e[i].b, e[i].x0, e[i].y0, e[i].pow);
  }

  if (!(GetArg(IniBuffer, "GenerateImage", Temp))) {
    GenerateImage = 1;
    Print(_DNormal, "'GenerateImage' entry missing in INI file. Setting it to %i\n", GenerateImage);
  } else {
    sscanf(Temp, "%i", &GenerateImage);
    Print(_DNormal, "'GenerateImage' entry is set to %i\n", GenerateImage);
  }

  if (!(GetArg(IniBuffer, "GenerateRadon", Temp))) {
    GenerateRadon = 1;
    Print(_DNormal, "'GenerateRadon' entry missing in INI file. Setting it to %i\n", GenerateRadon);
  } else {
    sscanf(Temp, "%i", &GenerateRadon);
    Print(_DNormal, "'GenerateRadon' entry is set to %i\n", GenerateRadon);
  }

  if (GenerateImage == 1) {
    ImImage = NewFloatImage("Image",
                            XSamples * (2 * OverSamp + 1),
                            YSamples * (2 * OverSamp + 1),
                            _RealArray);
    ZeroImage(ImImage);
  }

  if (GenerateRadon == 1) {
    RtImage = NewFloatImage("Radon",
                            ThetaSamples * (2 * OverSamp + 1),
                            RhoSamples * (2 * OverSamp + 1),
                            _RealArray);
    ZeroImage(RtImage);
  }
}

/**********************************************************
[NAME]
end

[SYNOPSIS]
void init(void)

[DESCRIPTION]
The function will write the images and release
the allocated memory.

[REVISION]
25/10-94 Peter Toft 
**********************************************************/
void end(void) {
  char t1[80], t2[80];
  int m, n, tooversampp1, i, j;
  float *poi1, *poi2;
  float sum;

  tooversampp1 = 2 * OverSamp + 1;
  strcpy(t1, OutFileName);
  strcpy(t2, OutFileName);
  strcat(t1, "a");
  strcat(t2, "r");

  if (GenerateImage == 1) {
    OutImImage = NewFloatImage("Image", XSamples, YSamples, _RealArray);
    OutImImage->DeltaX = DeltaX;
    OutImImage->DeltaY = DeltaX;
    OutImImage->Xmin = Xmin;
    OutImImage->Ymin = Ymin;
    RenameImage(OutImImage, t1);

    for (m = 0 ; m < XSamples ; m++) {
      poi2 = OutImImage->Signal[m];
      for (n = 0 ; n < YSamples ; n++) {
        sum = 0.0;
        for (i = m * tooversampp1 ; i < (m + 1) * tooversampp1 ; i++) {
          poi1 = ImImage->Signal[i];
          for (j = n * tooversampp1 ; j < (n + 1) * tooversampp1 ; j++) sum += poi1[j];
        }
        poi2[n] = sum / (tooversampp1 * tooversampp1);
      }
    }
    WriteImage(OutImImage, _FIF);
    FreeImage(OutImImage);
    FreeImage(ImImage);
  }

  if (GenerateRadon == 1) {
    OutRtImage = NewFloatImage("Radon", ThetaSamples, RhoSamples, _RealArray);
    OutRtImage->DeltaX = PI / ThetaSamples;
    OutRtImage->DeltaY = DeltaRho;
    OutRtImage->Xmin = 0.0;
    OutRtImage->Ymin = (1 - RhoSamples) * 0.5 * DeltaRho;
    RenameImage(OutRtImage, t2);

    for (m = 0 ; m < ThetaSamples ; m++) {
      poi2 = OutRtImage->Signal[m];
      for (n = 0 ; n < RhoSamples ; n++) {
        sum = 0.0;
        for (i = m * tooversampp1 ; i < (m + 1) * tooversampp1 ; i++) {
          poi1 = RtImage->Signal[i];
          for (j = n * tooversampp1 ; j < (n + 1) * tooversampp1 ; j++) sum += poi1[j];
        }
        poi2[n] = sum / (tooversampp1 * tooversampp1);
      }
    }
    WriteImage(OutRtImage, _FIF);
    FreeImage(RtImage);
    FreeImage(OutRtImage);
  }
  free(e);
}

/**********************************************************
[NAME]
imgen

[SYNOPSIS]
void imgen(void)

[DESCRIPTION]
The function will fill the image with the desired shapes.

[REVISION]
25/10-94 Peter Toft 
**********************************************************/
void imgen(void) {
  int i, xi, yi, LXSamples, LYSamples;
  double LDeltaX, LXmin, LYmin;
  float x, y;
  float *signal;

  LDeltaX = DeltaX / (2 * (float) OverSamp + 1.0);
  LXSamples = XSamples * (2 * OverSamp + 1);
  LYSamples = YSamples * (2 * OverSamp + 1);
  LXmin = Xmin - OverSamp * LDeltaX;
  LYmin = Ymin - OverSamp * LDeltaX;
  for (i = 0 ; i < NumberOfShapes ; i++) {
    Print(_DNoLog, "  shape number %3i (%3i)\r", i, NumberOfShapes);
    for (xi = 0 ; xi < LXSamples ; xi++) {
      x = xi * LDeltaX + LXmin;
      signal = ImImage->Signal[xi];
      for (yi = 0 ; yi < LYSamples ; yi++) {
        y = LYmin + yi * LDeltaX;
        signal[yi] += image_func(xi, yi, x, y, e[i]);
      }
    }
  }
  Print(_DNormal, "Signal has been generated\n");
  PrintStats(_DNormal, ImImage);
}

/**********************************************************
[NAME]
rtgen

[SYNOPSIS]
void rtgen(void)

[DESCRIPTION]
The function will fill the Radon space with the Radon
transform of the desired shapes.

[REVISION]
25/10-94 Peter Toft 
**********************************************************/
void rtgen(void) {
  int i, ti, ri;
  float theta, rholow, LDeltaTheta, LDeltaRho, thetalow;
  float *signal;
  float LThetaSamples, LRhoSamples;

  LThetaSamples = (2 * OverSamp + 1) * ThetaSamples;
  LRhoSamples = (2 * OverSamp + 1) * RhoSamples;
  LDeltaRho = DeltaRho / (2 * (float) OverSamp + 1.0);
  LDeltaTheta = PI / (float) LThetaSamples;
  rholow = (1 - RhoSamples) * 0.5 * DeltaRho - OverSamp * LDeltaRho;
  thetalow = -OverSamp * LDeltaTheta;
  for (i = 0 ; i < NumberOfShapes ; i++) {
    Print(_DNoLog, "  shape number %3i (%3i)\r", i, NumberOfShapes);
    for (ti = 0 ; ti < LThetaSamples ; ti++) {
      theta = ti * LDeltaTheta + thetalow;
      signal = RtImage->Signal[ti];
      for (ri = 0 ; ri < LRhoSamples ; ri++)
        signal[ri] += radon_func(rholow + ri * LDeltaRho, theta, e[i]);
    }
  }
  Print(_DNormal, "Radon domain has been generated\n");
  PrintStats(_DNormal, RtImage);
}

/**********************************************************
[NAME]
main

[SYNOPSIS]
void main(int argc,char *argv[])

[DESCRIPTION]
The main-function will control the program. 

[REVISION]
25/10-94 Peter Toft\\
April 12 96 Peter Toft
**********************************************************/
void main(int argc, char *argv[]) {
  idum = IntVector(1);
  idum[0] = -1;
  sqrt2 = sqrt(2.0);
  if (argc != 2) {
    Print(_DNormal, "Usage -> 'RadonAna inifile' to generate phantom\n");
    exit(0);
  }
  IniBuffer = ReadIni(argv[1]);      /* Read INI file into IniBuffer   */

  Print(_DNormal, "\n*******************************************\n\n");
  Print(_DNormal, "             RadonAna version 2.0\n");
  Print(_DNormal, "               by  Peter Toft \n");
  Print(_DNormal, "\n*******************************************\n");

  Print(_DNormal, "Initialize begin\n");
  init();
  Print(_DNormal, "Initialize done\n");
  if (GenerateImage == 1) {
    Print(_DNormal, "Image generation begin\n");
    imgen();
    Print(_DNormal, "Image generation done\n");
  }
  if (GenerateRadon == 1) {
    Print(_DNormal, "Radon generation begin\n");
    rtgen();
    Print(_DNormal, "Radon generation done\n");
  }
  Print(_DNormal, "Writing file\n");
  end();
}














