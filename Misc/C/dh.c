#include <stdio.h>
main()
{
    char c,d;
    while(!feof(stdin))
    {
	c=getc(stdin);
	if ((c==' ')||(c=='.')||(c==',')||(c=='?')||(c=='\n'))
	    putc(c,stdout);
	else
	{
	    d=getc(stdin);
	    if (d=='\n')
		putc(d,stdout);
	    else
	    {
		putc('~',stdout);
		putc('{',stdout);
		putc(c,stdout);
		putc(d,stdout);
		putc('~',stdout);
		putc('}',stdout);
	    }
	}
    }
    return;
}

