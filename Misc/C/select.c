/* cc -static select.c -o $BIN_PATH/select */

#include <stdio.h>
#include <string.h>
#define MAXCHAR_IN_LINE 500

int main(int argc,char **argv)
{
  char *head, *tail;
  char string[MAXCHAR_IN_LINE];
  char c;
  int i,j,BLANK_AS_HEAD=1,BLANK_AS_TAIL=1,NEW=1;
  
  if (argc>=2) 
   {
    BLANK_AS_HEAD=0;
    head = argv[1];
  }

  if (argc>=3) 
   {
    BLANK_AS_TAIL=0;
    tail = argv[2];
  }
  while(!feof(stdin))
   {
    i=readline(string);
    if ((i==0)&&BLANK_AS_HEAD) break;
    if ((!BLANK_AS_HEAD)&&(!strcmp(string,head))) break;
   }

  while(!feof(stdin))
   {
    i=readline(string);
    if ((!NEW)&&(i==0)&&BLANK_AS_TAIL) break;
    if ((!NEW)&&(!BLANK_AS_TAIL)&&(!strcmp(string,tail))) break;
    if (!(NEW&&(i=0))) 
      {
        NEW=0;
        puts(string);
      }
  }
return 0; 
}

int readline(char * string)
{
 int length=0;
 char *c=string;
 while (!feof(stdin))
   {
     if ((*c=getc(stdin))=='\n') break;
     else c++;
   }
 length = feof(stdin)?c-string-1:c-string;
 *c=0;
 return(length);
}

