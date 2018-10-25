/******************************************************************************
[HEADER] 

This library contains different common utilities used by most of the
other functions.

******************************************************************************/

#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <malloc.h>
#include "misc.h"
#include "imgtools.h"
#include "calc.h"

/*****************************************************************************
[NAME]
GetArg

[SYNOPSIS]
int GetArg(char *IniBuffer,
           char *Entry,
           char *Value);

[DESCRIPTION]

The function searches the text buffer {\tt IniBuffer} for the string 
{\tt 'Entry='}.  If the entry exists, the text on the right side of the
equal sign will be returned in {\tt Value} as a string.  If the
argument does not exist, the function will return with {\tt
\_ArgError}.

[USAGE]
{\tt Error=GetArg(IniBuffer, "Function", Function);}

Searches for an entry `Function' in {\tt IniBuffer}, and returns its value
in {\tt Function}.

[REVISION]
Sep. 94, JJJ and PAP
*******************************************************************************/ 
int GetArg(char *IniBuffer, char* Entry, char* Value)
{
  char *chrprt, *EntryPrt, *EntryEndPrt, SEntry[101], EEntry[101];
  int ArgError=_ArgError;

  SEntry[0]='\n';
  strcpy(&SEntry[1],Entry);
  strcpy(EEntry,SEntry);
  SEntry[(chrprt=strchr(SEntry,'\0'))-SEntry]=' ';
  chrprt[1]='\0';
  EEntry[(chrprt=strchr(EEntry,'\0'))-EEntry]='=';
  chrprt[1]='\0';
  if ( (chrprt=strstr(IniBuffer,SEntry)) || (chrprt=strstr(IniBuffer,EEntry)) ) {
    chrprt+=1; 
    if (!(EntryPrt=strchr(chrprt,'=')))      
      Error("Error in IniFile: '%s'", Entry);
    if (!(EntryEndPrt=strchr(chrprt,(char) 0x0A)))
      Error("Error in INI file: '%s'", Entry);
    strncpy(Value,(EntryPrt+=1),(EntryEndPrt-EntryPrt));
    Value[EntryEndPrt-EntryPrt]='\0';
    ArgError=_NoError; 
  }
  return(ArgError);
}


/*****************************************************************************
[NAME]
ReadIni

[SYNOPSIS]
char *ReadIni(char *FileName);

[DESCRIPTION]

This function finds and reads the file {\tt FileName.ini}, and copies
its contents into the return string. 

[USAGE]

\tc{ IniBuffer=ReadIni("Iradon.ini");}

Reads the {\tt .ini} file {\tt Iradon.ini} into {\tt IniBuffer}.

[REVISION]
Sep. 94, JJJ and PAP
*******************************************************************************/
char *ReadIni(char *FileName)
{
  char Temp[100], *chrprt, *IniBuf;
  FILE *InFile;
  int FileSize;

  strcpy(Temp,FileName);
  if(!(chrprt=strrchr(Temp,'.')))
    strcat(Temp,".ini");    

  if (!(InFile=fopen(Temp, "rb"))) 
    Error("Error opening INI file");

  fseek(InFile,0L,SEEK_END);
  FileSize=ftell(InFile);
  rewind(InFile);

  if (!(IniBuf=(char *)malloc(sizeof(char)*FileSize+3)))
    Error("Error allocating memory (ReadIni)");
  IniBuf[0]='\n';
  IniBuf[FileSize+1]='\n';
  IniBuf[FileSize+2]='\0'; 
  if (fread(&IniBuf[1], sizeof(char), FileSize, InFile)!=FileSize)
    Error("Error reading INI file"); 
  return IniBuf;
}


/*****************************************************************************
[NAME]
ReadIradonArgs

[SYNOPSIS]
void ReadIradonArgs(char *IniBuffer);

[DESCRIPTION]

This function searches the buffer {\tt IniBuffer} for all the standard
entries that are used to do the reconstruction with. The entries are
copied to the global structure {\tt IniFile}. Default values are used
if non-vital entries are not found.


[USAGE]

\tc{ReadIradonArgs(TestBuffer);}

Scans the buffer {\tt TestBuffer} for the standard entries.

[REVISION]
Sep. 94, JJJ and PAP\\
April 12, 96 PT typo
*****************************************************************************/
void ReadIradonArgs(char *IniBuffer)
{
  char Temp[100], *chrprt;
  Image *OrgImage;
 
  if (!(GetArg(IniBuffer,"OutFile",IniFile.OutFile)))
    Error("'OutFile' entry missing in INI file");
  chrprt=strrchr(IniFile.OutFile,'.');
  if(chrprt) {
    strcpy(Temp, chrprt+1);
    *chrprt='\0';
  }
       if (strstr(Temp,"gif")) IniFile.OutFileType=_GIF;
  else if (strstr(Temp,"fif")) IniFile.OutFileType=_FIF;
  else if (strstr(Temp,"dat")) IniFile.OutFileType=_DAT;
  else if (strstr(Temp,"mat")) IniFile.OutFileType=_MAT;
  else if (strstr(Temp,"analyze")) IniFile.OutFileType=_Analyze;
  else Error("Unknown file type: '%s'",Temp);

  if (!(GetArg(IniBuffer,"DebugLevel",IniFile.DebugLevel))) DebugNiveau=_DNormal;
  else if (strstr(IniFile.DebugLevel,"Normal"))   DebugNiveau=_DNormal;
  else if (strstr(IniFile.DebugLevel,"Debug"))    DebugNiveau=_DDebug;
  else if (strstr(IniFile.DebugLevel,"NoScreen")) DebugNiveau=_DNoScreen;
  else if (strstr(IniFile.DebugLevel,"NoLog"))    DebugNiveau=_DNoLog;
  else if (strstr(IniFile.DebugLevel,"HardCore")) DebugNiveau=_DHardCore;
  else Error("Unknown DebugLevel: '%s'",IniFile.DebugLevel);   

  OpenLog(IniFile.OutFile);

  if (!(GetArg(IniBuffer,"InFile",IniFile.InFile)))
    Error("'InFile' entry missing in INI file");
  chrprt=strrchr(IniFile.InFile,'.');
  if(chrprt) {
    strcpy(Temp, chrprt+1);
    *chrprt='\0';
  }
       if (strstr(Temp,"gif")) IniFile.InFileType=_GIF;
  else if (strstr(Temp,"fif")) IniFile.InFileType=_FIF;
  else if (strstr(Temp,"dat")) IniFile.InFileType=_DAT;
  else if (strstr(Temp,"mat")) IniFile.InFileType=_MAT;
  else if (strstr(Temp,"analyze")) IniFile.InFileType=_Analyze;
  else Error("Unknown file type: '%s'",Temp);

  if ((GetArg(IniBuffer,"OrgFile",IniFile.OrgFile))) { 
    chrprt=strrchr(IniFile.OrgFile,'.');
    if(chrprt) {
      strcpy(Temp, chrprt+1);
      *chrprt='\0';
    }
       if (strstr(Temp,"gif")) IniFile.OrgFileType=_GIF;
  else if (strstr(Temp,"fif")) IniFile.OrgFileType=_FIF;
  else if (strstr(Temp,"dat")) IniFile.OrgFileType=_DAT;
  else if (strstr(Temp,"mat")) IniFile.OrgFileType=_MAT;
  else if (strstr(Temp,"analyze")) IniFile.OrgFileType=_Analyze;
  else Error("Unknown file type: '%s'",Temp);
  }
  else strcpy(IniFile.OrgFile,"");

  if (!(GetArg(IniBuffer,"Function",IniFile.Function)))
    Error("'Function' entry missing in INI file");
 
  if (!(GetArg(IniBuffer,"Palette",IniFile.Palette)))
    strcpy(IniFile.Palette,"");

  if (!(GetArg(IniBuffer,"InterPol",Temp)))
    IniFile.InterPol=1;
  else
    sscanf(Temp,"%d",&(IniFile.InterPol));

  if (!(GetArg(IniBuffer,"FilterCutoff",Temp)))
    IniFile.FilterCutoff=1.0;
  else
    sscanf(Temp,"%f",&(IniFile.FilterCutoff));

  if (!(GetArg(IniBuffer,"FilterType",Temp)))
    IniFile.FilterType=_Ramp;
  else {
    sscanf(Temp,"%s",Temp);
    if (strstr(Temp,"Ramp")) IniFile.FilterType=_Ramp;
    else if (strstr(Temp,"Hanning")) IniFile.FilterType=_Hanning;
    else if (strstr(Temp,"Hamming")) IniFile.FilterType=_Hamming;
  }

   if (!(GetArg(IniBuffer,"SliceNumber",Temp)))
    IniFile.SliceNumber=1;
  else
    sscanf(Temp,"%d",&(IniFile.SliceNumber));
 
  if (strcmp(IniFile.OrgFile,"")) {
    OrgImage=ReadImage(IniFile.OrgFile,IniFile.OrgFileType);
    IniFile.Xmin=OrgImage->Xmin;
    IniFile.Ymin=OrgImage->Ymin;
    IniFile.XSamples=OrgImage->M;
    IniFile.YSamples=OrgImage->N;
    IniFile.DeltaX=OrgImage->DeltaX;
    IniFile.DeltaY=OrgImage->DeltaY;
    FreeImage(OrgImage);
  }
  else
    if ((strcmp(IniFile.Function,"Convert"))
	&& (strcmp(IniFile.Function,"Trace")) 
	&& (strcmp(IniFile.Function,"Info")) ){
       if (!(GetArg(IniBuffer,"Xmin",Temp)))
	sscanf(Temp,"%f",&(IniFile.Xmin));
      if (!(GetArg(IniBuffer,"Ymin",Temp)))
	Error("'Ymin' entry missing in INI file");
      sscanf(Temp,"%f",&(IniFile.Ymin));
      if (!(GetArg(IniBuffer,"Xmin",Temp)))
	Error("'Xmin' entry missing in INI file");
      sscanf(Temp,"%f",&(IniFile.Xmin));
      if (!(GetArg(IniBuffer,"DeltaX",Temp)))
	Error("'DeltaX' entry missing in INI file");
      sscanf(Temp,"%f",&(IniFile.DeltaX));
      if (!(GetArg(IniBuffer,"DeltaY",Temp)))
	Error("'DeltaY' entry missing in INI file");
      sscanf(Temp,"%f",&(IniFile.DeltaY));
      if (!(GetArg(IniBuffer,"XSamples",Temp)))
	Error("'XSamples' entry missing in INI file");
      sscanf(Temp,"%d",&(IniFile.XSamples));
      if (!(GetArg(IniBuffer,"YSamples",Temp)))
	Error("'YSamples' entry missing in INI file");
      sscanf(Temp,"%d",&(IniFile.YSamples));
  }
}


/*******************************************************************************
[NAME]
GetDateTime

[SYNOPSIS]
void GetDateTime(char *str,
                 int DateTimeFormat);
 
[DESCRIPTION]

This function gets the current time and date, and returns it in the
pointer {\tt str}. The return string are formatted accordingly to {\tt
DateTime} Format. Three formats are avaliable as of now, {\tt
\_LongDate}, {\tt \_Time} and \tc{RealTime}. {\tt \_Longdate} returns the date and
time in a string as `Day, dd. month, hh:mm:ss', {\tt \_Time} will
return `hh:mm:ss' and \tc{ \_RealTime} returns the number of seconds
from last midnight.

[USAGE]
{\tt GetDateTime(str,\_Time);}

Returns the current time in {\tt str} formatted as ``12:34:56''.

[REVISION]
Oct. 94, JJJ
*******************************************************************************/
void GetDateTime(char *str, int DateTimeFormat)
{
  time_t *timep;
  struct tm *times;

  if (!(timep=(time_t *)malloc(sizeof(time_t))))
    Error(" Memory allocation error (GetDateTime)");
/*  if (!(times=(struct tm *)malloc(sizeof(struct tm))))
    Error(" Memory allocation error (GetDateTime)");*/
  time(timep);
  times=localtime(timep);
  if (DateTimeFormat==_LongDate) 
    strftime(str, 55, "%A, %d. %B, %H.%M:%S", times);
  if (DateTimeFormat==_Time) 
    strftime(str, 20,"%H.%M:%S", times);
  if (DateTimeFormat==_RealTime)
    sprintf(str,"%d",times->tm_sec+times->tm_min*60+times->tm_hour*3600);
  Free(timep);
}


/*******************************************************************************
[NAME]
OpenLog

[SYNOPSIS]
void OpenLog(char* FileName);

[DESCRIPTION]

This function creates and opens a log file named {\tt
FileName.log}. Appends time and date to the top if the file.  If the
global variable {\tt DebugNiveau} equals {\tt \_DHardCore} or {\tt
\_DNolog} no log will be opened. The file handle is the global
variable {\tt LogFile}.

[USAGE]
{\tt OpenLog("Test.log");}

Creates and opens the log file {\tt Test.log}.

[REVISION]
Oct. 94, JJJ
*******************************************************************************/
void OpenLog(char *FileName)
{
  char str[100];
  if (!(DebugNiveau&(_DHardCore|_DNoLog))) {  
    strcpy(LogFileName,FileName);
    strcat(LogFileName,".log");
    GetDateTime(str, _LongDate);
    if (!(LogFile=fopen(LogFileName, "w"))) 
      Error("Error opening LOG file");
    fprintf(LogFile,"--- START OF IRADON LOG; %s ---\n", str);   
    fclose(LogFile);
  }
}


/*******************************************************************************
[NAME]
CloseLog

[SYNOPSIS]
void CloseLog()

[DESCRIPTION]

This function closes the log file and appends time and date to the end
of the file.

[USAGE]
{\tt CloseLog():}

Closes the current log file.

[REVISION]
Oct. 94, JJJ
*******************************************************************************/

void CloseLog()
{
  char str[100];

  if (!(DebugNiveau&(_DHardCore|_DNoLog))) {  
    fclose(LogFile);
    GetDateTime(str, _Time);
    if (!(LogFile=fopen(LogFileName, "a+"))) 
      Error("Error opening LOG file");
    fprintf(LogFile,"---- Normal termination; %s ---\n", str);   
    fclose(LogFile);
  }
}


/*******************************************************************************
[NAME]
Print

[SYNOPSIS]
void Print(int Niveau,
           char *fmt,
           ...);

[DESCRIPTION]

This function handles all the writing to the screen and to the
log. {\tt fmt} and {\tt ...} are the same arguments as to {\tt
print}. The variable {\tt \_Niveau} combined with the global variable
{\tt DebugNiveau} chooses what should be written. {\tt DebugNiveau}
takes precedence over {\tt Niveau}. The following choises exist

\begin{tabbing}
{\tt \_DHardCore} \== Nothing written.\\
{\tt \_DSilent}   \>== Nothing written on screen, full log.\\
{\tt \_DNolog}    \>= Full information on screen, no log.\\
{\tt \_DNormal}   \>= Some information on screen, full log.\\
{\tt \_DDebug}    \>= Full information both on screen and in log.
\end{tabbing}

[USAGE]
\verb+Print(_DNormal,"File not found '%%s' \n", TestFileName);+ 

Prints the above message both on the screen and in the log file.

[REVISION]
Oct. 94, JJJ
*******************************************************************************/
void Print(int Niveau, char *fmt, ...)
{ 
  char LogString[255];
  char str[20];
  va_list ap;

  va_start(ap,fmt);
  vsprintf(LogString, fmt, ap);  
  GetDateTime(str, _Time);
  if ((Niveau!=_DNoLog) && (!(DebugNiveau&(_DHardCore|_DNoLog)))) {
    if (!(LogFile=fopen(LogFileName, "a+"))) 
      Error("Error opening LOG file");
    fprintf(LogFile,"%s: %s",str,LogString);
    fclose(LogFile);
  }
  if (((DebugNiveau&(_DDebug|_DNoLog)) && (Niveau&(_DDebug|_DNoLog|_DNormal))) 
     || ((DebugNiveau&(_DNormal)) && (Niveau&(_DNormal|_DNoLog)))) {
    vprintf(fmt,ap);
    fflush(stdout);
  }
  va_end(ap);
} 


/*******************************************************************************
[NAME]
Error

[SYNOPSIS]
void Error(char *fmt,
           ...);

[DESCRIPTION]

This function is called when a fatal error occurs. It prints the error
message on the screen, closes the log and exits. {\tt fmt} and {\tt
...} are the same arguments as to {\tt print}.

[USAGE]
{\tt Error("Memory allocation problems");}

Prints the above message, closes the log and exits.


[REVISION]
Oct. 94, JJJ and PAP
*******************************************************************************/
void Error(char *fmt, ...)
{
  char LogString[255];
  char str[100];
  va_list ap;

  va_start(ap,fmt);
  vsprintf(LogString, fmt, ap); 
  va_end(ap);

  if (!(DebugNiveau&(_DHardCore|_DNoLog))) {
    GetDateTime(str, _Time);
    if (!(LogFile=fopen(LogFileName, "a+"))) 
      printf("Error opening LOG file");
    fprintf(LogFile,"%s: %s\n",str,LogString);
    fprintf(LogFile,"---- IRADON terminated with an error; %s ---\n", str);   
    fclose(LogFile);
  }

  printf(LogString);
  printf("\n");
  exit(1);
}


/**********************************************************
[NAME]
MultReStore

[SYNOPSIS]
void MultReStore(float *p1,float *p2)

[DESCRIPTION]
The function will multiply the two complex numbers $A$
and $B$ and store the result at the location of $A$. Here
$A=p1[0]~+~i~p1[1]$ and $B=p2[0]~+~i~p2[1]$.

[USAGE]

{\t tMultReStore(arr1,arr2);}

Preforms the multiplication {\tt arr1=arr1*arr2}.

[REVISION]
Nov. 94, JJJ and PT
**********************************************************/
void MultReStore(float *p1,float *p2)
{
  multtemp=p1[0]*p2[0]-p1[1]*p2[1];
  p1[1]=p1[0]*p2[1]+p1[1]*p2[0];
  p1[0]=multtemp;
}


/**********************************************************
[NAME]
MultNew

[SYNOPSIS]
void MultNew(float *p1,float *p2,float *p3)

[DESCRIPTION] The function will multiply the two complex numbers $A$
and $B$, so that $C=AB$. The result is stored at the location of
$C$. Here $A=p1[0]~+~i~p1[1]$, $B=p2[0]~+~i~p2[1]$ and
$C=p3[0]~+~i~p3[1]$.

[USAGE]
{\tt MultNew(arr1,arr2,arr3);}

Preforms the multiplication {\tt arr3=arr1*arr2}.

[REVISION]
Nov. 94, JJJ and PT 
**********************************************************/
void MultNew(float *p1,float *p2,float *p3)
{
  p3[0]=p1[0]*p2[0]-p1[1]*p2[1];
  p3[1]=p1[0]*p2[1]+p1[1]*p2[0];
}


/********************************************************************************
[NAME]
FloatVector

[SYNOPSIS]
float *FloatVector(int Size);

[DESCRIPTION]

Allocates and returns a vector of floats, with \tc{Size} number of
elements. The vector is initialized to 0. Checking is made to ensure
no memory problems.

[USEAGE]
{\tt Test=FloatVector(99);}

Allocates \tc{Test} as a float vector 99 elements long.

[REVISION]
Oct. 94, JJJ
*****************************************************************************/
float *FloatVector(int Size)
{
  float *data;

  if (!(data=(float*)calloc(Size,sizeof(float))))
    Error("Memory allocation problems (FloatVector). %i elements.",Size);
  return data;
}


/*****************************************************************************
[NAME]
IntVector

[SYNOPSIS]
float *IntVector(int Size);

[DESCRIPTION]

Allocates and returns a vector of integers, with \tc{Size} number of
elements. The vector is initialized to 0. Adequate checking is done.

[USEAGE]
{\tt Test=IntVector(99);}

Allocates \tc{Test} as an int vector 99 elements long.

[REVISION]
Oct. 94, JJJ
*****************************************************************************/
int *IntVector(int Size)
{
  int *data;

  if (!(data=(int*)calloc(Size,sizeof(int))))
    Error("Memory allocation problems (IntVector). %i elements",Size);
  return data;
}


/*****************************************************************************
[NAME]
convert2

[SYNOPSIS]
void convert2(char *ind);

[DESCRIPTION]

The function swaps the two bytes at the location {\tt ind}. Used for
converting short int's between HP/PC.

[USEAGE]
{\tt convert2((char *)\&X);}

Converts {\tt X}, where {\tt X} is a short int.

[REVISION]
Oct. 94, PAP and PT
*****************************************************************************/
void convert2(char *ind)
{
  char tmp;
 
  tmp=ind[0];
  ind[0]=ind[1];
  ind[1]=tmp;
}


/******************************************************************************
[NAME]
convert4

[SYNOPSIS]
void convert4(char *ind);

[DESCRIPTION]

The function swaps the four bytes at the location {\tt ind} around the
center. Used for converting short floats between HP/PC.

[USEAGE]
{\tt convert4((char *)\&X);}

Converts {\tt X}, where {\tt X} is a float.

[REVISION]
Oct. 94, PAP and PT
*****************************************************************************/
void convert4(char *ind)
{
  char tmp;

  tmp=ind[0];
  ind[0]=ind[3];
  ind[3]=tmp;
  tmp=ind[1];
  ind[1]=ind[2];
  ind[2]=tmp;
}


/*****************************************************************************
[NAME]
convert8

[SYNOPSIS]

void convert8(char *ind);

[DESCRIPTION]

The function swaps the eight bytes at the location {\tt ind} around
the center. Used for converting doubles between HP/PC.

[USEAGE]
{\tt convert8((char *)\&X);}

Converts {\tt X}, where {\tt X} is a double.

[REVISION]
Oct. 94, PAP and PT
*****************************************************************************/
void convert8(char *ind)
{
  char tmp;

  tmp=ind[0];
  ind[0]=ind[7];
  ind[7]=tmp;
  tmp=ind[1];
  ind[1]=ind[6];
  ind[6]=tmp;
  tmp=ind[2];
  ind[3]=ind[5];
  ind[5]=tmp;
  tmp=ind[3];
  ind[3]=ind[4];
  ind[4]=tmp;
}


/*****************************************************************************
[NAME]
convert\_UNIX\_PC

[SYNOPSIS]
void convert_UNIX_PC(char *ind, int lgd);

[DESCRIPTION]

The function swaps the {\tt lgd} bytes at the location {\tt ind}
around the center. Used for converting arbitrary length numbers between
HP/PC.

[USEAGE]
{\tt convert\_UNIX\_PC((char *)\&X, 8);}

Converts {\tt X}, where {\tt X} is 8 bytes long eg. a double.

[REVISION]
Oct. 94, PAP and PT
****************************************************************************/
void convert_UNIX_PC(char *ind, int lgd)
{
  int i,lgdh,lgdm1;
  char tmp;

  if (lgd>1)
  {
    lgdh=lgd/2;
    lgdm1=lgd-1;
    for (i=0;i<lgdh;i++)
    {
      tmp=ind[i];
      ind[i]=ind[lgdm1-i];
      ind[lgdm1-i]=tmp;
    }
  }
}


/*****************************************************************************
[NAME]
strequal

[SYNOPSIS]
int strequal(char *str1, char *str2);

[DESCRIPTION]
This function will compare two strings by length and see whether one is
contained in the other.

[USEAGE]
{\tt ok=strequal("Radon","Radon transform");}
This will return a zero. Drop transform and it will return a 1.

[REVISION]
April 9, 96 PToft
****************************************************************************/
int strequal(char *str1, char *str2)
{
  if ((strlen(str1)==strlen(str2)) && (strstr(str1,str2)==str1))
    return 1;
  else
    return 0;
}
