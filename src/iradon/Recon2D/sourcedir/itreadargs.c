/* itgetpar.c */

#include "it.h"

char Algorithmar[3][50] = {"Algebraric Reconstruction",
                           "Expection Maximization",
                           "Conjugate Gradients"};

char RadonKernelar[6][50] = {"Sinc interpolation",
                             "Binary Nearest neighboar",
                             "Ray driven nearest neighboar",
                             "Ray driven linear interpolation",
                             "Pixel driven (rays)",
                             "Pixel driven (for all points)"};

char TrueFalsear[2][6] = {"False",
                          "True"};

extern itINItype itINI;

/*****************************************************************************
[NAME]
ReadItArgs

[SYNOPSIS]
void ReadIt(char *IniBuffer);

[DESCRIPTION]

This function searches the buffer {\tt IniBuffer} for all the standard
entries that are used to do the reconstruction with. The entries are
copied to the global structure {\tt IniFile}. Default values are used
if non-vital entries are not found.


[USAGE]

\tc{ReadItArgs(TestBuffer);}

Scans the buffer {\tt TestBuffer} for the standard entries.

[REVISION]
April 96 JJJ and PT
*****************************************************************************/
void ReadItArgs(char *IniBuffer) {
  char Temp[100];
  char *strptr;
  Image *tempImage, *refImage = NULL;

  if (!(GetArg(IniBuffer, "OutFileName", Temp)))
    Error("`OutFileName' entry missing in INI file");
  sscanf(Temp, "%s", itINI.OutFileName);
  if (!(strptr = strrchr(itINI.OutFileName, '.')))
    Error("Error in IniFile: `%s'", itINI.OutFileName);
  itINI.OutFileName[strptr - itINI.OutFileName] = '\0';
  Print(_DNormal, "`OutFileName' entry is set to `%s'\n", itINI.OutFileName);

  if (!(GetArg(IniBuffer, "DebugLevel", Temp))) DebugNiveau = _DNormal;
  else if (strstr(Temp, "Normal")) DebugNiveau = _DNormal;
  else if (strstr(Temp, "Debug")) DebugNiveau = _DDebug;
  else if (strstr(Temp, "NoScreen")) DebugNiveau = _DNoScreen;
  else if (strstr(Temp, "NoLog")) DebugNiveau = _DNoLog;
  else if (strstr(Temp, "HardCore")) DebugNiveau = _DHardCore;
  else Error("Unknown DebugLevel: `%s'", Temp);

  strcpy(LogFileName, itINI.OutFileName);
  OpenLog(LogFileName);

  if (!(GetArg(IniBuffer, "InFileName", Temp)))
    Error("`InFileName' entry missing in INI file");
  sscanf(Temp, "%s", itINI.InFileName);
  if (!(strptr = strrchr(itINI.InFileName, '.')))
    Error("Error in IniFile: `%s'", itINI.InFileName);
  itINI.InFileName[strptr - itINI.InFileName] = '\0';
  Print(_DNormal, "`InFileName' entry is set to `%s'\n", itINI.InFileName);

  tempImage = ReadFIFHeader(itINI.InFileName);
  Print(_DNormal, "`ThetaSamples' entry is set to %i\n", (itINI.ThetaSamples = tempImage->M));
  Print(_DNormal, "`ThetaMin' entry is set to %f\n", (itINI.ThetaMin = tempImage->Xmin));
  Print(_DNormal, "`DeltaTheta' entry is set to %f\n", (itINI.DeltaTheta = tempImage->DeltaX));
  Print(_DNormal, "`RhoSamples' entry is set to %i\n", (itINI.RhoSamples = tempImage->N));
  Print(_DNormal, "`RhoMin' entry is set to %f\n", (itINI.RhoMin = tempImage->Ymin));
  Print(_DNormal, "`DeltaRho' entry is set to %f\n", (itINI.DeltaRho = tempImage->DeltaY));
  Free(tempImage);

  if ((GetArg(IniBuffer, "RefFileName", Temp))) {
    sscanf(Temp, "%s", itINI.RefFileName);
    if (!(strptr = strrchr(itINI.RefFileName, '.')))
      Error("Error in IniFile: `%s'", itINI.RefFileName);
    itINI.RefFileName[strptr - itINI.RefFileName] = '\0';
    Print(_DNormal, "RefImage = '%s.fif'\n", itINI.RefFileName);
    refImage = ReadFIFHeader(itINI.RefFileName);
    Print(_DNormal, "Output sampling parameters taken from refImage: %s.fif'\n",
          itINI.RefFileName);
    Print(_DNormal, "`XSamples' entry is set to %i\n", (itINI.XSamples = refImage->M));
    Print(_DNormal, "`YSamples' entry is set to %i\n", (itINI.YSamples = refImage->N));
    Print(_DNormal, "`Xmin' entry is set to %f\n", (itINI.Xmin = refImage->Xmin));
    Print(_DNormal, "`Ymin' entry is set to %f\n", (itINI.Ymin = refImage->Ymin));
    Print(_DNormal, "`DeltaX' entry is set to %f\n", (itINI.DeltaX = refImage->DeltaX));
    Print(_DNormal, "`DeltaY' entry is set to %f\n", (itINI.DeltaY = itINI.DeltaX));
    Free(refImage);
  } else {
    Print(_DNormal, "`RefFileName' entry missing in INI file.  Feature not used\n");
    strcpy(itINI.RefFileName, "\0");
  }

  if ((GetArg(IniBuffer, "StartFileName", Temp))) {
    sscanf(Temp, "%s", itINI.StartFileName);
    if (!(strptr = strrchr(itINI.StartFileName, '.')))
      Error("Error in IniFile: `%s'", itINI.StartFileName);
    itINI.StartFileName[strptr - itINI.StartFileName] = '\0';
    Print(_DNormal, "Startguess = '%s.fif'\n",
          itINI.StartFileName);
    tempImage = ReadFIFHeader(itINI.StartFileName);
    if (strlen(itINI.RefFileName) != 0)
      if (memcmp(&refImage->M, &tempImage->M, sizeof(float) * 7) != 0)
        Error("RefImage and StartImage have different sampling parameters.");
      else Print(_DNormal, "ref- and startImage has equal sampling parameters.\n");
    if (strlen(itINI.RefFileName) == 0) {
      Print(_DNormal, "Output sampling parameters taken from startImage, `%s'\n",
            itINI.StartFileName);
      Print(_DNormal, "`XSamples' entry is set to %i\n", (itINI.XSamples = tempImage->M));
      Print(_DNormal, "`YSamples' entry is set to %i\n", (itINI.YSamples = tempImage->N));
      Print(_DNormal, "`Xmin' entry is set to %f\n", (itINI.Xmin = tempImage->Xmin));
      Print(_DNormal, "`Ymin' entry is set to %f\n", (itINI.Ymin = tempImage->Ymin));
      Print(_DNormal, "`DeltaX' entry is set to %f\n", (itINI.DeltaX = tempImage->DeltaX));
      Print(_DNormal, "`DeltaY' entry is set to %f\n", (itINI.DeltaY = itINI.DeltaX));
    }
    if (strlen(itINI.RefFileName) != 0)
      Free(refImage);
    Free(tempImage);
  } else {
    Print(_DNormal, "`StartFileName' entry missing in INI file.  Feature not used\n");
    strcpy(itINI.StartFileName, "\0");
  }
  if (strlen(itINI.RefFileName) == 0 && strlen(itINI.StartFileName) == 0) {
    Print(_DNormal, "Reading output sampling parameters from `INI' file.\n");
    if (!(GetArg(IniBuffer, "XSamples", Temp)))
      Error("`XSamples' entry missing in INI file");
    sscanf(Temp, "%i", &itINI.XSamples);
    Print(_DNormal, "`XSamples' entry is set to %i\n", itINI.XSamples);

    if (!(GetArg(IniBuffer, "YSamples", Temp)))
      Print(_DNormal, "`YSamples' entry missing in INI file. Setting it to XSamples=%i\n",
            (itINI.YSamples = itINI.XSamples));
    else {
      sscanf(Temp, "%i", &itINI.YSamples);
      Print(_DNormal, "`YSamples' entry is set to %i\n", itINI.YSamples);
    }

    if (!(GetArg(IniBuffer, "DeltaX", Temp)))
      Error("`DeltaX' entry missing in INI file");
    sscanf(Temp, "%f", &itINI.DeltaX);
    Print(_DNormal, "`DeltaX' entry is set to %f\n", itINI.DeltaX);
    itINI.DeltaY = itINI.DeltaX;
    Print(_DNormal, "`DeltaY' equals 'DeltaX', set to %f\n", itINI.DeltaY);

    if (!(GetArg(IniBuffer, "Xmin", Temp)))
      Print(_DNormal, "`Xmin' entry missing in INI file. Setting it to %f\n",
            (itINI.Xmin = (1 - itINI.XSamples) * 0.5 * itINI.DeltaX));
    else {
      sscanf(Temp, "%f", &itINI.Xmin);
      Print(_DNormal, "`Xmin' entry is set to %f\n", itINI.Xmin);
    }

    if (!(GetArg(IniBuffer, "Ymin", Temp)))
      Print(_DNormal, "`Ymin' entry missing in INI file. Setting it to %f\n",
            (itINI.Ymin = (1 - itINI.YSamples) * 0.5 * itINI.DeltaX));
    else {
      sscanf(Temp, "%f", &itINI.Ymin);
      Print(_DNormal, "`Ymin' entry is set to %f\n", itINI.Ymin);
    }
  }
  if (!(GetArg(IniBuffer, "LowestALevel", Temp)))
    Print(_DNormal, "`LowestALevel' entry missing in INI file. Setting it to %f\n",
          (itINI.LowestALevel = 0.0));
  else {
    sscanf(Temp, "%f", &itINI.LowestALevel);
    Print(_DNormal, "`LowestALevel' entry is set to %f\n", itINI.LowestALevel);
  }

  if (!(GetArg(IniBuffer, "ConstrainMin", Temp)))
    Print(_DNormal, "`ConstrainMin' entry missing in INI file. Setting it to %f\n",
          (itINI.ConstrainMin = -1.0));
  else {
    sscanf(Temp, "%f", &itINI.ConstrainMin);
    Print(_DNormal, "`ConstainMin' entry is set to %f\n", itINI.ConstrainMin);
  }

  if (!(GetArg(IniBuffer, "ConstrainMax", Temp)))
    Print(_DNormal, "`ConstrainMax' entry missing in INI file. Setting it to %f\n",
          (itINI.ConstrainMax = -1.0));
  else {
    sscanf(Temp, "%f", &itINI.ConstrainMax);
    Print(_DNormal, "`ConstainMax' entry is set to %f\n", itINI.ConstrainMax);
  }

  if (!(GetArg(IniBuffer, "Iterations", Temp)))
    Error("`Iterations' entry missing in INI file");
  sscanf(Temp, "%i", &itINI.Iterations);
  Print(_DNormal, "`Iterations' entry is set to %i\n", itINI.Iterations);

  if (!(GetArg(IniBuffer, "SaveIterations", Temp)))
    Print(_DNormal, "`SaveIterations' entry missing in INI file. Setting it to %i\n",
          (itINI.SaveIterations = 0));
  else {
    sscanf(Temp, "%i", &itINI.SaveIterations);
    if (itINI.SaveIterations) {
      if (itINI.Iterations / itINI.SaveIterations > MaxSaves) {
        Print(_DNormal, "To large number of saves requested, aborting SaveIterations\n");
        itINI.SaveIterations = 0;
      } else
        Print(_DNormal, "Saving temporary result every %i iterations\n",
              itINI.SaveIterations);
    } else Print(_DNormal, "No temporary saving \n");
  }
  if (!(GetArg(IniBuffer, "Alpha", Temp)))
    Print(_DNormal, "`Alpha' entry missing in INI file. Setting it to %f\n",
          (itINI.Alpha = 1.0));
  else {
    sscanf(Temp, "%f", &itINI.Alpha);
    Print(_DNormal, "`Alpha' entry is set to %f\n", itINI.Alpha);
  }

  if (!(GetArg(IniBuffer, "Beta", Temp)))
    Print(_DNormal, "`Beta' entry missing in INI file. Setting it to %f\n",
          (itINI.Beta = 1.0));
  else {
    sscanf(Temp, "%f", &itINI.Beta);
    Print(_DNormal, "`Beta' entry is set to %f\n", itINI.Beta);
  }

  if (!(GetArg(IniBuffer, "Regularization", Temp)))
    Print(_DNormal,
          "`Regularization' entry missing in INI file. Setting it to %f\n",
          itINI.Regularization = 0);
  else {
    sscanf(Temp, "%f", &itINI.Regularization);
    Print(_DNormal,
          "`Regularization' entry is set to %f\n",
          itINI.Regularization);
  }

  if (!(GetArg(IniBuffer, "KernelFileName", Temp)))
    Error("`KernelFileName' entry missing in INI file");
  sscanf(Temp, "%s", itINI.KernelFileName);
  if (!(strptr = strstr(itINI.KernelFileName, ".sif")))
    Error("Error in IniFile: `%s'", itINI.KernelFileName);
  itINI.KernelFileName[strptr - itINI.KernelFileName] = '\0';
  Print(_DNormal, "`KernelFileName' entry is set to `%s'\n", itINI.KernelFileName);

  if (!(GetArg(IniBuffer, "SaveMatLab", Temp)))
    Print(_DNormal, "`SaveMatLab' entry missing in INI file. Setting it to %s\n", TrueFalsear[(itINI.SaveMatLab = 0)]);
  else {
    sscanf(Temp, "%d", &itINI.SaveMatLab);
    Print(_DNormal, "`SaveMatLab' entry is set to %s\n", TrueFalsear[itINI.SaveMatLab]);
  }

  if (!(GetArg(IniBuffer, "IterationType", Temp)))
    Error("`IterationType' entry missing in INI file");
  sscanf(Temp, "%i", &itINI.IterationType);
  Print(_DNormal, "`Numbering' entry is set to %i\n", itINI.IterationType);

  if (!(GetArg(IniBuffer, "Algorithm", Temp))) Error("`Algorithm' must be specified.");
  else if (strcmp(Temp, "EM") == 0) itINI.Algorithm = _EM;
  else if (strcmp(Temp, "ART") == 0) itINI.Algorithm = _ART;
  else if (strcmp(Temp, "CG") == 0) itINI.Algorithm = _CG;
  else Error("Unknown Algorithm: `%s'", Temp);
  Print(_DNormal, "`Algorithm' entry is set to %s\n", Algorithmar[itINI.Algorithm]);

  if (!(GetArg(IniBuffer, "UseFast", Temp))) {
    Print(_DNormal, "No speed setting, using FAST method");
    itINI.IsFast = 1;
  }
  sscanf(Temp, "%i", &itINI.IsFast);
  Print(_DNormal, "`UseFast' entry is set to %i\n", itINI.IsFast);

  if (!(GetArg(IniBuffer, "RadonKernel", Temp)))
    Error("`RadonKernel' must be specified.");

  else if (strcmp(Temp, "SINC") == 0) itINI.RadonKernel = _SINC;
  else if (strcmp(Temp, "NN") == 0) itINI.RadonKernel = _NN;
  else if (strcmp(Temp, "RNN") == 0) itINI.RadonKernel = _RNN;
  else if (strcmp(Temp, "RL") == 0) itINI.RadonKernel = _RL;
  else if (strcmp(Temp, "P1") == 0) itINI.RadonKernel = _P1;
  else if (strcmp(Temp, "P2") == 0) itINI.RadonKernel = _P2;
  else Error("Unknown RadonKernel: `%s'", Temp);
  Print(_DNormal, "`RadonKernel' entry is set to %s\n",
        RadonKernelar[itINI.RadonKernel]);

  if (!(GetArg(IniBuffer, "OverSamp", Temp)))
    Print(_DNormal, "`OverSamp' entry missing in INI file. Setting it to %s\n", itINI.OverSamp = 0);
  else {
    sscanf(Temp, "%d", &itINI.OverSamp);
    Print(_DNormal, "`OverSamp' entry is set to %d\n", itINI.OverSamp);
  }
}









