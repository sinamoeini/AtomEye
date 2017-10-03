#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define buf_size 1024*1024 /* 1 M */ 

main (int argc,char **argv)
{
 FILE * in;
 FILE * out;
 char *buffer;
 char *end;
 char *begin;
 char *p;
 char format[15];
 int i,first=1,last,Para_finish=1;
 long size=buf_size;
 if (argc==2)
 {
  first=1;
  last=atoi(argv[1]);
  }
 else
 {
  first=atoi(argv[1]);
  last=atoi(argv[2]);
  }
 if (last==first)
  {
   sprintf(format,"uudecode %d ; rm %d",first,first);
   system(format);
   return;
   } 
 while (!(buffer=malloc(sizeof(char)*size))) size/=2;
 out=fopen ("uncoded.xrn","w+");
 for (i=first;i<=last;i++)
  {
   sprintf(buffer,"%d",i);
   if (!(in=fopen(buffer,"r")))
    {
     printf("Unable to open file %s...exit.\n",buffer);
     exit(0);
     }
   while(1)
    {
     begin=buffer; 
     fread(buffer,size,1,in);
     end=buffer+ftell(in);
     if (Para_finish&&(i==first))
      {
       while(1)
        {
         if (!(begin=strstr(begin,"begin")))
          {
           printf("No begin line in the file %d...exit.\n",i);
           exit(1);
           }
         if (*(begin-1)=='\n') break;
         begin++; 
         }
        }
      else 
       if (Para_finish)
        while(1)
         {
          while(*begin!='M') begin++;
          if ((*(begin+61)=='\n')&&(*(begin+62)='M')) break;
          begin++;
	  }
       if ((!Para_finish)&&(i==last)) begin=buffer;
       p=begin;
       if (Para_finish) while(*(p++)!='\n');
        else while (*p!='M') p++;
       while((*p=='M')&&(p<end)) p+=62;
       if ((p>=end)||(i==last)) p=end;
       if (fwrite(begin,p-begin,1,out)!=1)
        {
         printf("Not enough disk space...exit.\n");
         exit(2);
         }
       if(Para_finish=feof(in)) break;
       }
      fclose(in);
      }
     fclose(out);
     free(buffer);
     if (last<10) 
      {
       sprintf(format,"rm [%d-%d]",first,last);
       system(format);
       }
     else 
       {
        sprintf(format,"rm [%d-9]",first);
        system(format);
        for(i=1;i<last/10;i++)
	 {
          sprintf(format,"rm %d[0-9]",i);
          system(format);
	  }
        sprintf(format,"rm %d[0-%d]",i,last-last/10);
        system(format);
        }        
     system("uudecode uncoded.xrn");
     return;
   }  




































