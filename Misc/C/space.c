#include <stdio.h>

main()
{
 char c,d;
 while(1)
   {
    c=getc(stdin);
    if (feof(stdin)) break;
    if ( c=='\n') 
       { 
         putc('\n',stdout);
         d = getc(stdin);
         if (d!='\n') 
	   {
	     putc('\n',stdout);
	     putc(d,stdout);
	   }
       }
    else putc(c,stdout);
  }

return;
}
