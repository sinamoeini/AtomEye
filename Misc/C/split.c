#include <stdio.h>
#include <stdlib.h>

int main(int argc,char **argv)
{
 long bytes=atol(argv[2]);
 char s[10];
 int a=1;
 FILE * m;
 FILE * v;
 long c;
 void * w = malloc(bytes*sizeof(char));
 v=fopen(argv[1],"r");
 while(1)
  {
   sprintf (s,"%d",a);
   m=fopen(s,"w");
   c=fread(w,1,bytes,v);
   fwrite(w,c,1,m);
   fclose(m);
   printf("%s\n",s);
   if (c!=bytes) break;
   a++;
   }
 fclose(v);
 return 1;
 }





