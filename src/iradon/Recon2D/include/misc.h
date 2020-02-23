/* Misc.h */

typedef struct {
  char InFile[100];
  int InFileType;
  char OrgFile[100];
  int OrgFileType;
  char OutFile[100];
  int OutFileType;
  char Function[50];
  char DebugLevel[50];
  float Param[10];
  char Palette[100];
  float Xmin;
  float Ymin;
  float DeltaX;
  float DeltaY;
  int XSamples;
  int YSamples;
  int InterPol;
  int FilterType;
  float FilterCutoff;
  int SliceNumber;
} INI;

#define _DHardCore  1
#define _DNoScreen  2
#define _DNoLog     4
#define _DNormal    8
#define _DDebug    16

#define _LongDate 100
#define _Time     101
#define _RealTime 102

#define max(a, b) ((a)>(b) ? (a) : (b))
#define min(a, b) ((a)<(b) ? (a) : (b))
#define sq(a) ((a)*(a))

extern int DebugNiveau;
extern FILE *LogFile;
extern char LogFileName[100];
extern INI IniFile;
extern int DebugNiveau;
extern char LogFileName[100];
extern float multtemp;

extern int GetArg(char *IniBuffer, char *Entry, char *Value);
extern char *ReadIni(char *FileName);
extern void ReadIradonArgs(char *IniBuffer);
extern void GetDateTime(char *str, int DateTimeFormat);

extern void OpenLog(char *);
extern void CloseLog(void);
extern void Print(int, char *, ...);
extern void Error(char *, ...);

extern void MultNew(float *p1, float *p2, float *p3);
extern void MultReStore(float *p1, float *p2);
extern float *FloatVector(int Size);
extern int *IntVector(int Size);

extern void convert2(char *);
extern void convert4(char *);
extern void convert8(char *);
extern void convert_UNIX_PC(char *ind, int lgd);
int strequal(char *, char *);










