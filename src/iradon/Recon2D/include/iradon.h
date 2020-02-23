#include <math.h>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <malloc.h>
#include "imgtools.h"
#include "calc.h"
#include "eval.h"
#include "misc.h"

int DebugNiveau = _DNoLog;
FILE *LogFile;
char LogFileName[100];
float multtemp;
INI IniFile;


