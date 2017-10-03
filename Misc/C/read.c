#include <stdio.h>
#include <stdlib.h>

main()
{
  char c;
  char name[200] = "uudecode        ";
  char rm[200] = "rm ";
  char *p = name + 9;
  char *q = rm + 3;

  while (!feof(stdin))
    {
      c = getc(stdin);
      if ((c!=' ')&&(c!='\n'))
	{
	  if ((c>='0')&&(c<='9'))
	    {
	      while ((c!='\t')&&(c!='\n'))
		{
		  *(q++)=*(p++) = c;
		  c = getc(stdin);
		}
	      *p = 0;
	      system (name);
	      system (rm);
	      p = name+9;
	      q = rm+3;
	    }
	  else
	    while ((getc(stdin)!=' ') && !feof(stdin));
	}
    }
}
