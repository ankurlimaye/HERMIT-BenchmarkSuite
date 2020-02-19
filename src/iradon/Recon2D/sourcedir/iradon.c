/***************************************************************************
[HEADER]

This program is the main program, which contains routines to call the
different direct reconstruction routines. The mode of operation
in determined by the {\tt Function=} entry in the {\tt .ini} file the
program is called with.

All functions in this module are not intended for outside use, they are
called accordingly to the contents of the used {\tt .ini} file.

The ordering of the following items are not important, but each require
on line.
\begin{description}
\item[{\tt OutFile}] Name of the file to write the reconstructed image (for
the reconstruction functions). The
program recognizes the following file extentions: 
\begin{description}
\item[{\tt fif}] Float image format which includes sampling parameters. 
\item[{\tt gif}] The very popular graphics format. 
\item[{\tt dat}] raw and outdated imageformat. 
\item[{\tt mat}] Matlab files with one matrix. 
\item[{\tt analyze}] GE specific file format.
\end{description}
\item[{\tt InFile}] Name of the file to read the sinogram image (for
the reconstruction functions). The sinogram should be a fif-file and
contain the sampling parameters.
\item[{\tt OrgFile}] (Optional) If this is provided as a fif-file the
sampling parameters are used to create the output-file, and measures
of misfit can be made between the OrgFile and the OutFile.
\item[{\tt Function}] The function the program should do. Currently
supported functions are
\begin{description}
\item[{\tt FB}] Filtered Backprojection.
\item[{\tt BF}] Filtering After Backprojection.
\item[{\tt CNF}] Central Slice. FFT based with Nearest Neighbor approximation.
\item[{\tt CBF}] Central Slice. FFT based with Bilinear Interpolation.
\item[{\tt CNC}] Central Slice. Chirp-z based with Nearest Neighbor approximation.
\item[{\tt CBC}] Central Slice. Chirp-z based with Bilinear Interpolation.
\item[{\tt Forward}] The InFile is forward Radon transformed into the OutFile.
\item[{\tt Convert}] Very handy function Convert Infile directly to OutFile
(in another file-format).
\item[{\tt Trace}] Pick a trace of the InFile.
\item[{\tt Info}] Get statistics about the Infile.
\end{description}
\item[{\tt DebugLevel}] This parameter controls the level of output. The
parameter is mixed and overrules with the one used in the
Print-statements.
\begin{description}
\item[{\tt Normal}] Standard level of output to screen and logfile.
\item[{\tt Debug}] Allmost all output is logged to screen and logfile.
\item[{\tt NoScreen}] Screen output is disabled and logfile level is ordinary.
\item[{\tt NoLog}] Logfile output is disabled and screen level is ordinary.
\item[{\tt HardCore}] No information at all.
\end{description}
\item[{\tt Palette}] (Optional) Name of external colormap file to use when writing gif-images.
\item[{\tt InterPol}] (Optional) Interpolation level. Used by Filtered Backprojection.
\item[{\tt SliceNumber}] (Optional) If the InFile is an analyze-file with a lot of
slices, then this parameter will select the slice number.
\item[{\tt Xmin}] The minimum x-position of the reconstructed image. Not
needed if a fif-file is provided as OrgFile.
\item[{\tt Ymin}] The minimum y-position of the reconstructed image.  Not
needed if a fif-file is provided as OrgFile.
\item[{\tt DeltaX}] Sampling distance on the x-axis.  Not needed if a
fif-file is provided as OrgFile.
\item[{\tt DeltaY}] Sampling distance on the y-axis. Not needed if a fif-file
is provided as OrgFile.
\item[{\tt XSamples}] Number of samples on the x-axis. Not needed if a
fif-file is provided as OrgFile.
\item[{\tt YSamples}] Number of samples on the y-axis. Not needed if a
fif-file is provided as OrgFile.
\end{description}

Assume the the sinogram {\tt mysino.fif} should be reconstructed into
{\tt rec.fif} using
Filtered Backprojection with the same parameters as found in {\tt
myorg.fif} then the parameter file {\tt myrec.ini} could look like:
{\small\begin{verbatim}
InFile=mysino.fif
OutFile=rec.fif
Function=FB
DebugLevel=Normal
InterPol=1
OrgFile=myorg.fif
\end{verbatim}}


PT April 96
***************************************************************************/

#include "iradon.h"
char *IniBuffer;

void DoForward(void) {
  Image *MyImage, *MySino;

  MyImage = ReadImage(IniFile.InFile, IniFile.InFileType);
  MySino = ReadImage(IniFile.OrgFile, IniFile.OrgFileType);

  Forward(MyImage, MySino);
  RenameImage(MySino, IniFile.OutFile);
  WriteImage(MySino, IniFile.OutFileType);
  FreeImage(MyImage);
  FreeImage(MySino);
}

/***************************************************************************

[NAME]
DoFilteredBack

[SYNOPSIS]
DoFilteredBack();

[DESCRIPTION] 

This function reads the image specified in the \tc{.ini}-file, performes an
reconstruction using Filtered Backprojection, renames the image and
saves the result. Calls \tc{FilteredBack}.

[USAGE] 
{\tt DoFilteredBack();} 

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void DoFilteredBack(void) {
  Image *NewImage, *InvNewImage, *ref;

  Print(_DNormal, "Using Filtered Backprojection\n");
  NewImage = ReadImage(IniFile.InFile, IniFile.InFileType);
  InvNewImage = FilteredBack(NewImage);
  if (strcmp(IniFile.OrgFile, "")) {
    ref = ReadImage(IniFile.OrgFile, IniFile.OrgFileType);
    Print(_DNormal, "L2 = %9.6f \n", L2Norm(ref, InvNewImage));
    FreeImage(ref);
  }
  RenameImage(InvNewImage, IniFile.OutFile);
  WriteImage(InvNewImage, IniFile.OutFileType);
  FreeImage(NewImage);
  FreeImage(InvNewImage);
}

/***************************************************************************
[NAME]
DoBackFiltering

[SYNOPSIS]
DoBackFiltering();

[DESCRIPTION] 

This function reads the image specified in the \tc{.ini}-file,
performes reconstruction using Filtering after Backprojection, renames
the image and saves the result. Calls \tc{BackFiltering}.

[USAGE] 

{\tt DoBackFiltering();} 

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void DoBackFiltering(void) {
  Image *NewImage, *InvNewImage, *ref;

  Print(_DNormal, "Using Filtering after Backprojection\n");
  NewImage = ReadImage(IniFile.InFile, IniFile.InFileType);
  InvNewImage = BackFilter(NewImage);
  if (strcmp(IniFile.OrgFile, "")) {
    ref = ReadImage(IniFile.OrgFile, IniFile.OrgFileType);
    Print(_DNormal, "L2 = %9.6f \n", L2Norm(ref, InvNewImage));
    FreeImage(ref);
  }
  RenameImage(InvNewImage, IniFile.OutFile);
  WriteImage(InvNewImage, IniFile.OutFileType);
  FreeImage(NewImage);
  FreeImage(InvNewImage);
}

/***************************************************************************
[NAME]
DoCentralSliceChirp

[SYNOPSIS]
DoCentralSliceChirp();

[DESCRIPTION] 

This function reads the image specified in the \tc{ini}-file, performs
reconstruction using Central Slice with either nearest neighbor or
chirp-z interpolation, renames the image and saves the result. Calls
\tc{CentralSliceNN}, \tc{CentralSliceBL} or
\tc{CentralSliceCZ}. Inverts the spectrum using inverse two
dimensional chirp-z transformation.

[USAGE] 

{\tt DoCentralSliceChirp();} 

[REVISION]
Oct. 94, JJJ and PAP
***************************************************************************/
void DoCentralSliceChirp(void) {
  char Value[100];
  Image *NewImage, *spectrum, *ref;

  Print(_DNormal, "Using Central Slice with Chirp-z interpolation\n");
  NewImage = ReadImage(IniFile.InFile, IniFile.InFileType);
  GetArg(IniBuffer, "Function", Value);
  if (strequal(Value, "CNC"))
    spectrum = CentralSliceNN(NewImage);
  else if (strequal(Value, "CBC"))
    spectrum = CentralSliceBL(NewImage);
  else
    spectrum = CentralSliceCZ(NewImage);
  FreeImage(NewImage);

  NewImage = IChirpSpectrum(spectrum);
  if (strcmp(IniFile.OrgFile, "")) {
    ref = ReadImage(IniFile.OrgFile, IniFile.OrgFileType);
    Print(_DNormal, "L2 = %9.6f \n", L2Norm(ref, NewImage));
    FreeImage(ref);
  }
  RenameImage(NewImage, IniFile.OutFile);
  WriteImage(NewImage, IniFile.OutFileType);
  FreeImage(NewImage);
  FreeImage(spectrum);
}

/***************************************************************************
[NAME]
DoCentralSliceFFT

[SYNOPSIS]
DoCentralSliceFFT();

[DESCRIPTION] 

This function read the image specified in the \tc{ini}-file, performes
reconstruction using Central Slice with either nearest neighbor or
bi-linear interpolation, renames the image and saves the
result. Calls \tc{CentralSliceNN}, \tc{CentralSliceBL} or
\tc{CentralSliceCZ}. Uses the two dimensional FFT to invert the
spectrum.

[USAGE] 

{\tt DoCentralSliceChirp();} 

[REVISION]
Oct. 94, JJJ and PAP\\
April 96 PToft More info is written and better check on Value
***************************************************************************/
void DoCentralSliceFFT(void) {
  char Value[100];
  Image *NewImage, *spectrum = NULL, *ref;

  NewImage = ReadImage(IniFile.InFile, IniFile.InFileType);
  GetArg(IniBuffer, "Function", Value);
  if (strequal(Value, "CNF")) {
    Print(_DNormal,
          "Using Central Slice with FFT - Nearest Neighbor interpolation\n");
    spectrum = CentralSliceNN(NewImage);
  } else if (strequal(Value, "CBF")) {
    Print(_DNormal,
          "Using Central Slice with FFT - Bilinear interpolation\n");
    spectrum = CentralSliceBL(NewImage);
  } else
    Error("Bad function (DoCentralSliceFFT)\n");

  NewImage = IFFTSpectrum(spectrum);
  if (strcmp(IniFile.OrgFile, "")) {
    ref = ReadImage(IniFile.OrgFile, IniFile.OrgFileType);
    Print(_DNormal, "L2 = %9.6f \n", L2Norm(ref, NewImage));
    FreeImage(ref);
  }
  RenameImage(NewImage, IniFile.OutFile);
  WriteImage(NewImage, IniFile.OutFileType);
  FreeImage(NewImage);
  FreeImage(spectrum);
}

/***************************************************************************
[NAME]
DoConvert

[SYNOPSIS]
DoConvert();

[DESCRIPTION] 

This function reads the image specified in the \tc{ini}-file, renames
the image and saves in the format specified by the {\tt Outfile} entry.

[USAGE] 

{\tt DoConvert();} 

[REVISION]
Nov. 94, JJJ and PAP
***************************************************************************/
void DoConvert(void) {
  Image *NewImage;
  char temp1[100], temp2[100];
  float min, max;
  int ABS;

  NewImage = ReadImage(IniFile.InFile, IniFile.InFileType);
  RenameImage(NewImage, IniFile.OutFile);

  if ((IniFile.InFileType == _GIF) && (IniFile.OutFileType == _GIF))
    WriteGif(NewImage, "", 8, 0, 255, FALSE);
  else {
    if (((GetArg(IniBuffer, "GifMin", temp1) == _NoError)
        && (GetArg(IniBuffer, "GifMax", temp2) == _NoError))) {
      if ((temp1[0] != '\0') && (temp2[0] != '\0')) {
        sscanf(temp1, "%f", &min);
        sscanf(temp2, "%f", &max);
        if (min > max)
          Error("Invalid `GifMin' or `GifMax' specified");
        if ((GetArg(IniBuffer, "GifAbs", temp1) == _NoError)) {
          if (temp1[0] != '\0') {
            sscanf(temp1, "%i", &ABS);
            if ((ABS != 0) && (ABS != 1))
              ABS = 0;
          }
        } else ABS = 0;
      } else {
        min = 0.0;
        max = 0.0;
      }
    }

    if (IniFile.OutFileType == _GIF)
      WriteGif(NewImage, "", 8, min, max, ABS);
    else
      WriteImage(NewImage, IniFile.OutFileType);
  }
  FreeImage(NewImage);
}

/***************************************************************************
[NAME]
DoInfo

[SYNOPSIS]
DoInfo();

[DESCRIPTION] 

This function reads the image specified in the \tc{ini}-file, and
prints some information about it, including size and sampling
parameters. Calls the function {\tt PrintStats}.

[USAGE]

{\tt DoInfo();} 

[REVISION]
Nov. 94, JJJ
***************************************************************************/
void DoInfo(void) {
  Image *NewImage;

  NewImage = ReadImage(IniFile.InFile, IniFile.InFileType);
  Print(_DNoLog, "\n");
  PrintStats(_DNoLog, NewImage);
  Print(_DNoLog, "\n");
  FreeImage(NewImage);
}

/***************************************************************************
[NAME]
DoTrace

[SYNOPSIS]
DoTrace();

[DESCRIPTION] 

This function reads the image specified in the \tc{ini}-file, and
writes an ASCII-file containing the values in the picture for a given
row or collumn. The file is readable with gnuplot. Uses the entries
{\tt LineNr} and {\tt Dir} in the {\tt .ini} file to pick the row or
collumn.

[USAGE]

{\tt DoTrace();} 

[REVISION]
Nov. 94, JJJ
***************************************************************************/
void DoTrace(void) {
  int LineNr, Dir;
  char DirStr[100], TraceName[100], Temp[100];
  Image *NewImage;
  NewImage = ReadImage(IniFile.InFile, IniFile.InFileType);

  GetArg(IniBuffer, "Dir", Temp);
  sscanf(Temp, "%s", DirStr);
  Dir = ((strchr(DirStr, 'v')) ? _Vertical : _Horizontal);

  if (GetArg(IniBuffer, "LineNr", Temp)) {
    sscanf(Temp, "%d", &LineNr);
    if ((LineNr < 0) || (LineNr > ((Dir == _Horizontal) ? NewImage->M : NewImage->N)))
      LineNr = NewImage->N / 2;
  }

  if (GetArg(IniBuffer, "TraceFile", Temp))
    sscanf(Temp, "%s", TraceName);
  else
    strcpy(TraceName, IniFile.OutFile);

  printf("Tracing: %s %s %d \n", TraceName, ((Dir == _Horizontal) ? "Horz" : "Vert"), LineNr);
  WriteTrace(NewImage, TraceName, LineNr, Dir);
  FreeImage(NewImage);
}

/***************************************************************************
[NAME]
DoTestChirpz

[SYNOPSIS]
DoTestChirpz();

[DESCRIPTION] 

Internal debugging routine.

[USAGE]

{\tt DoTestChirpz();} 

[REVISION]
Nov. 94, JJJ
***************************************************************************/
void DoTestChirpz(void) {
  Image *i;
  int n;

  i = NewFloatImage("vector", 1, 128, _ComplexArray);

  for (n = 48 ; n < 80 ; n++)
    i->Signal[0][2 * n] = 1.0;

  WriteTrace(i, "signal", 0, _Vertical);
  ComplexFFT(i->Signal[0], i->N, _FFT);
  WriteTrace(i, "signal.1", 0, _Vertical);
  ComplexFFT(i->Signal[0], i->N, _IFFT);
  WriteTrace(i, "signal.2", 0, _Vertical);

  FreeImage(i);
  i = NewFloatImage("vector", 1, 128, _ComplexArray);
  for (n = 48 ; n < 80 ; n++)
    i->Signal[0][2 * n] = 1.0;
  ComplexFFT(i->Signal[0], i->N, _FFT);
  FFTShift(i, 1);
  i->Signal[0] = ComplexChirpZ(i->Signal[0], i->N, i->N,
                               1.0 / (i->N) * 34, 0.2 / (i->N), _IFFT);
  for (n = 0 ; n < 128 ; n++) {
    i->Signal[0][2 * n] /= 128.0;
    i->Signal[0][2 * n + 1] /= 128.0;
  }
  RealImage(i);
  WriteTrace(i, "signal.3", 0, _Vertical);
  FreeImage(i);
}

/***************************************************************************
[NAME]
DoTest

[SYNOPSIS]
DoTest();

[DESCRIPTION] 

Function to run all the avaliable reconstruction routines on the same
image. The $L_1$ and $L_2$ measures are calculated for each routine.

[USAGE]

{\tt DoTest();} 

[REVISION]
Nov. 94, JJJ and PAP
***************************************************************************/
void DoTest(void) {
  Image *Image1, *Image2, *ImageFB, *ImageBF, *ImageCSBL,
      *ImageCSNN, *ImageCSCH, *spectrum, *ref;
  float Tid;

  Image1 = ReadImage(IniFile.InFile, IniFile.InFileType);
  Image2 = CopyImage(Image1);
  ref = ReadImage(IniFile.OrgFile, IniFile.OrgFileType);

  Print(_DNormal, "1: \n");
  Tid = clock();
  spectrum = CentralSliceNN(Image1);
  ImageCSNN = IFFTSpectrum(spectrum);
  Tid = (clock() - Tid) / (float) CLOCKS_PER_SEC;
  Print(_DNormal, "IRadon: Reconstruction was active for %.2f seconds\n", Tid);
  FreeImage(Image1);
  FreeImage(spectrum);

  Print(_DNormal, "2:                    \n");
  Image1 = CopyImage(Image2);
  Tid = clock();
  spectrum = CentralSliceBL(Image1);
  ImageCSBL = IFFTSpectrum(spectrum);
  Tid = (clock() - Tid) / (float) CLOCKS_PER_SEC;
  Print(_DNormal, "IRadon: Reconstruction was active for %.2f seconds\n", Tid);
  FreeImage(Image1);
  FreeImage(spectrum);

  Print(_DNormal, "3:                    \n");
  Image1 = CopyImage(Image2);
  Tid = clock();
  spectrum = CentralSliceCZ(Image1);
  ImageCSCH = IChirpSpectrum(spectrum);
  Tid = (clock() - Tid) / (float) CLOCKS_PER_SEC;
  Print(_DNormal, "IRadon: Reconstruction was active for %.2f seconds\n", Tid);
  FreeImage(Image1);
  FreeImage(spectrum);

  Print(_DNormal, "4:                    \n");
  Image1 = CopyImage(Image2);
  Tid = clock();
  ImageBF = BackFilter(Image1);
  NormImage(ImageBF, 1.0, MeanValue(ref));
  Tid = (clock() - Tid) / (float) CLOCKS_PER_SEC;
  Print(_DNormal, "IRadon: Reconstruction was active for %.2f seconds\n", Tid);
  FreeImage(Image1);

  Print(_DNormal, "5:                    \n");
  Image1 = CopyImage(Image2);
  Tid = clock();
  ImageFB = FilteredBack(Image1);
  Tid = (clock() - Tid) / (float) CLOCKS_PER_SEC;
  Print(_DNormal, "IRadon: Reconstruction was active for %.2f seconds\n", Tid);
  FreeImage(Image1);
  FreeImage(Image2);

  Print(_DNormal, "FB:   L1=%9.6f, L2=%9.6f \n", L1Norm(ref, ImageFB), L2Norm(ref, ImageFB));
  Print(_DNormal, "BF:   L1=%9.6f, L2=%9.6f \n", L1Norm(ref, ImageBF), L2Norm(ref, ImageBF));
  Print(_DNormal, "CSNN: L1=%9.6f, L2=%9.6f \n", L1Norm(ref, ImageCSNN),
        L2Norm(ref, ImageCSNN));
  Print(_DNormal, "CSBL: L1=%9.6f, L2=%9.6f \n", L1Norm(ref, ImageCSBL),
        L2Norm(ref, ImageCSBL));
  Print(_DNormal, "CSCH: L1=%9.6f, L2=%9.6f \n", L1Norm(ref, ImageCSCH),
        L2Norm(ref, ImageCSCH));
/*
  Print(_DNormal,"Writing images...\n");
  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,"_cb");
  RenameImage(ImageCSBL,TempName);
  WriteImage(Diff,IniFile.OutFileType);
  Diff=DiffImage(ref,ImageCSBL);
  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".diff.csbl");
  RenameImage(Diff,TempName);
  WriteImage(Diff,IniFile.OutFileType);
  FreeImage(Diff);
  
  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".image.csnn");
  RenameImage(ImageCSNN,TempName);
  Diff=DiffImage(ref,ImageCSNN);
  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".diff.csnn");
  RenameImage(Diff,TempName);
  WriteImage(Diff,IniFile.OutFileType);
  FreeImage(Diff);

  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".image.csch");
  RenameImage(ImageCSCH,TempName);
  Diff=DiffImage(ref,ImageCSCH);
  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".diff.csch");
  RenameImage(Diff,TempName);
  WriteImage(Diff,IniFile.OutFileType);
  FreeImage(Diff);

  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".image.bf");
  RenameImage(ImageBF,TempName);
  Diff=DiffImage(ref,ImageBF);
  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".diff.bf");
  RenameImage(Diff,TempName);
  WriteGif(Diff,"",8,-1.0,1.0,TRUE);  
  WriteImage(Diff,IniFile.OutFileType);
  FreeImage(Diff);

  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".image.fb");
  RenameImage(ImageFB,TempName);
  Diff=DiffImage(ref,ImageFB);
  strcpy(TempName,IniFile.OutFile);
  strcat(TempName,".diff.fb");
  RenameImage(Diff,TempName);
  WriteGif(Diff,"",8,-1.0,1.0,TRUE);  
  WriteImage(Diff,IniFile.OutFileType);
  FreeImage(Diff);
*/

  FreeImage(ImageFB);
  FreeImage(ImageBF);
  FreeImage(ImageCSNN);
  FreeImage(ImageCSBL);
  FreeImage(ImageCSCH);

  FreeImage(ref);
}

/***************************************************************************
[NAME]
Main

[SYNOPSIS]
ProgrameName IniFile.ini

[DESCRIPTION] 

This is the main function. It reads the specified {\tt IniFile.ini}
file and executes one of the following functions:

\begin{itemize}
\item CentralSlice.
\item FilteredBack.
\item BackFiltering.
\item Convert.
\item Trace.
\item ImageInfo.
\item Test.
\end{itemize}

[USAGE] 

{\tt ProgramName Test.ini} 

Executes the program and uses {\tt Test.ini} as input.

[REVISION]
Oct. 94, JJJ and PAP\\
April 5, 96 PT More info\\
April 9, 96 PT Better check on Value\\
April 13, 96 PT Error in percent print if 0 sec.
***************************************************************************/
int main(int argc, char *argv[]) {
  int RealTid1, RealTid2;
  float Tid;
  char Value[100];

  GetDateTime(Value, _RealTime);
  sscanf(Value, "%i", &RealTid1);
  Tid = clock();

  if (!(argc == 2))
    Error("iradon, : You must specify ini-file");
  IniBuffer = ReadIni(argv[1]);
  ReadIradonArgs(IniBuffer);

  Print(_DNormal, "---------------------------\n");
  Print(_DNormal, "      iradon (ver 2.0)     \n");
  Print(_DNormal, "                           \n");
  Print(_DNormal, " Made by Jesper J. Jensen  \n");
  Print(_DNormal, "         Peter Philipsen   \n");
  Print(_DNormal, "         Peter Toft        \n");
  Print(_DNormal, "                           \n");
  Print(_DNormal, "---------------------------\n");

  GetArg(IniBuffer, "Function", Value);
  if (strequal(Value, "CC") || strequal(Value, "CNC") || strequal(Value, "CBC"))
    DoCentralSliceChirp();
  else if (strequal(Value, "CNF") || strequal(Value, "CBF"))
    DoCentralSliceFFT();
  else if (strequal(Value, "FB"))
    DoFilteredBack();
  else if (strequal(Value, "BF"))
    DoBackFiltering();
  else if (strequal(Value, "Forward"))
    DoForward();
  else if (strequal(Value, "Convert"))
    DoConvert();
  else if (strequal(Value, "Info"))
    DoInfo();
  else if (strequal(Value, "Trace"))
    DoTrace();
  else if (strequal(Value, "Test"))
    DoTest();
  else
    Error("Function not recognized: '%s'", Value);

  GetDateTime(Value, _RealTime);
  sscanf(Value, "%i", &RealTid2);
  Tid = (clock() - Tid) / (float) CLOCKS_PER_SEC;
  Print(_DNormal, "iradon: Program was active for %.2f seconds\n", Tid);
  Print(_DNormal, "        World time elapsed %d seconds\n",
        (RealTid2 - RealTid1));
  if (RealTid2 != RealTid1)
    Print(_DNormal, "        Program used %.2f %% cpu time\n",
          Tid / ((float) RealTid2 - RealTid1) * 100);
  Print(_DNormal, "-----------------------\n");
  CloseLog();
  return 0;
}
