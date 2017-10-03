#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

main (int argc,char **argv)
{
 int first,last;
 char format[20];
 int second = 8;
 
 if (argc==2) second = atoi(argv[1]);

 while(1)
   {
     sleep(second);	 
     system ("ls | /mit/liju99/Src/Bin/read");
   }
 return;
}

































