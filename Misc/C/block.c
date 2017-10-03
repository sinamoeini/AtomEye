#include <stdio.h>
#include <string.h>
#define MAXCHAR 100000

int main(int argc, char** argv)
{
 char *a,*b,*c;
 char buffer[MAXCHAR];
 fread (buffer,MAXCHAR,1,stdin);
 a = strstr(buffer,argv[1]);
 b = strstr(a,argv[2]);
 c = strchr(b, '\n');
 fwrite (a,c-a+1,1,stdout);
 return(0);
}   
