#include <stdio.h>

main()
{
 char c;
 while(!feof(stdin))
   {
    c=getc(stdin);
    if ( (c<='z')&&(c>='a')) putc(c+'A'-'a',stdout);
     else putc(c,stdout);
  }
return;
}
