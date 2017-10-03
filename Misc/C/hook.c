#include <stdio.h>

main(int argc,char **argv)
{
 char s[300];
 int hook=1;
 int turn=atoi(argv[1]);
 long j=1;
 if (argc==3) hook=atoi(argv[2]);
 while(!feof(stdin))
  {
   fgets(s,300,stdin);
   if (j%turn==hook) fputs(s,stdout);
   j++;
   }
  return;
}

