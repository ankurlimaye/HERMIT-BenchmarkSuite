#include <stdio.h>
#include <string.h>

FILE *infilep,*outfilep;
char *poi,*EndFile,*chrprt,outfilename[200],*infilename,text[200],cfilename[200];
int hea,syn,des,not,rev,nam,usa,inc,incomment,ii,jj;

int main(int argc,char * argv[])
{
  printf("\nC comment to LaTeX translator is running\n\n");
  
  if (argc!=2)
  {
    printf("This program requires one argument! Use CLAT CFILE (- .c)\n");
    exit(1);
  }
  infilename=(char *)malloc(200);
  strcpy(infilename,argv[1]);
  strcat(infilename,".c");
  infilep=fopen(infilename,"rt");
  if (infilep==NULL)
  {
    printf("The file %s does not exist\n",infilename);
    exit(1);
  }
  if (strrchr(infilename,'/')) 
  {
    chrprt=strrchr(infilename,'/');
    chrprt+=1;
  }
  else
    chrprt=infilename;

  jj=0;
  for (ii=0;ii<strlen(chrprt);ii++)
  {
    if (chrprt[ii]=='_')
      text[jj++]='\\';
    text[jj++]=chrprt[ii];
  }
  text[jj]='\0';

  strcpy(outfilename,argv[1]);
  strcat(outfilename,".tex");
  outfilep=fopen(outfilename,"wt");
  fprintf(outfilep,"\\section{%s}\n",text);
  EndFile=fgets(text,200,infilep);
  incomment=0;
  
  while (EndFile!=NULL)
  {
    if ((strstr(text,"/**")==text) || (incomment==1))
    {
      incomment=1;
      if (strstr(text,"[")==text)
      {
        nam=(strstr(text,"[NAM")==text) ? 1 : 0;
        syn=(strstr(text,"[SYN")==text) ? 1 : 0;
        des=(strstr(text,"[DES")==text) ? 1 : 0;
        usa=(strstr(text,"[USA")==text) ? 1 : 0;
        not=(strstr(text,"[NOT")==text) ? 1 : 0;
        rev=(strstr(text,"[REV")==text) ? 1 : 0;
        inc=(strstr(text,"[INC")==text) ? 1 : 0;
        hea=(strstr(text,"[HEA")==text) ? 1 : 0;

        if (hea==1) {
	  fgets(text,200,infilep);
	  while ((strstr(text,"[")!=text) && (strstr(text,"**/")==NULL))
	  {
	    fprintf(outfilep,text);
	    fgets(text,200,infilep);
	  }
        }        
        else 
	if (syn==1)
	{
          fprintf(outfilep,"\n\\item[Synopsis]\n\\hspace*{2cm}\n");
          fprintf(outfilep,"\\vspace{\\itemsep}\\vspace{-\\baselineskip}\n\\begin{verbatim}\n");
	  fgets(text,200,infilep);
	  while ((strstr(text,"[")!=text) && (strstr(text,"**/")==NULL))
	  {
	    if (strstr(text,"\n")!=text)
	      fprintf(outfilep,text);
	    fgets(text,200,infilep);
    	  }	  
          fprintf(outfilep,"\\end{verbatim}\n\\vspace{\\itemsep}\\vspace{-\\baselineskip}\n");
	}
	else
	{
	  if (nam==1)
	  {
	    fgets(text,200,infilep);
	    printf("Found function %s",text);
            poi=strrchr(text,'\n');
	    if (poi!=NULL) *poi='\0';          /* change last return til slash 0 */
       	    fprintf(outfilep,"\n\\subsection{%s}\n",text);
       	    fprintf(outfilep,"\\hspace*{1em}\n");
	    fprintf(outfilep,"\\begin{minipage}[t]{16cm}\n");
	    fprintf(outfilep,"\\begin{description}");
	  }
	  else
	  {
	    if (des==1) fprintf(outfilep,"\\item[Description]\n\\hspace*{1mm}\n\n");
	    if (not==1) fprintf(outfilep,"\n\\item[Note]\n\\hspace*{1mm}\n\n");
	    if (inc==1) fprintf(outfilep,"\n\\item[Includes]\n\\hspace*{1mm}\n\n");
	    if (usa==1) fprintf(outfilep,"\n\\item[Usage]\n\\hspace*{1mm}\n\n");
	    if (rev==1) fprintf(outfilep,"\n\\item[Revision]\n\\hspace*{1mm}\n\n");
	    EndFile=fgets(text,200,infilep);
	    while ((strstr(text,"[")!=text) && (strstr(text,"**/")==NULL))
	    {
	      fprintf(outfilep,text);
	      fgets(text,200,infilep);
	    }	  
	  }
	}
      }
      else
	fgets(text,200,infilep);
      if ((strstr(text,"**/")!=NULL)&&(hea!=1))
      {
	incomment=0;
	fprintf(outfilep,"\\end{description}\n\\end{minipage}\n");
      }
    }
    else
      EndFile=fgets(text,200,infilep);
  }
  
  fclose(infilep);
  fclose(outfilep);
  printf("\nThe translator has written the LaTeX file %s\n\n",outfilename);
  return 0;
}













