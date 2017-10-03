/*
  Copyright (C) 1989      Fung F. Lee

  sgb2hz: convert a Macintosh/CCDOS SGB file into a HZ file.

  This program is free for general distribution.  

  This program runs on UNIX. You are welcome to port it to other operating
  systems.

*/

#include <stdio.h>

#define true 1
#define false 0
#define notAscii(c)	((c)&0x80)

int CR2LF=false;	/* flag for converting ASCII <CR> to <LF> */
int termStyle = false;	/* enforce terminal emulation style if true */
int MAXLEN = 77;	/* default maximum line length in the above style */
int MINLEN = 7;		/* minimum line length in the above style */
char *progname;

main(argc, argv)
     int argc;
     char *argv[];
{
  int i, len;

  progname = argv[0];
  for (i=1; i<argc; i++)
    {
      if (argv[i][0]!='-') warning();
      switch (argv[i][1])
	{
	case 'n': CR2LF = true; break;
	case 't': termStyle = true;
	  if ((argv[i][2]!='\0') && (sscanf(&argv[i][2], "%d", &len) != 1))
	    warning();
	  if (len >= MINLEN) MAXLEN = len;
	  break;
	default:  warning(); break;
	}
    }
  filter(stdin, stdout);
}

warning()
{
  fprintf(stderr, "usage: %s [-n] [-t[maximum-line-length]] < foo > foo.hz\n",
	  progname);
  exit(1);
}


filter(fin, fout)
     FILE *fin, *fout;
{
  int c1, c2, c3, c4, hi;
  int GBmode = false;
  int len = 0;
  
  while ((c1=fgetc(fin)) != EOF)
    {
      if (notAscii(c1))
#ifdef MAC
	{
	  hi = c1 & 0xF0;
	  switch (hi)
	    {
	    case 0x80:
	    case 0x90:
	    case 0xA0:
	      if (termStyle)
		{
		  if (GBmode && len>MAXLEN-5)
		    {
		      fprintf(fout, "~}~\n");
		      GBmode = false; len = 0;
		    }
		  else if (!GBmode && len>MAXLEN-7)
		    {
		      fprintf(fout, "~\n");
		      GBmode = false; len = 0;
		    }
		}
	      if (!GBmode) /* switch to GB mode */
		{
		  fprintf(fout, "~{");
		  len += 2;
		}
	      GBmode = true;
	      c2 = fgetc(fin);
	      mac2gb(c1, c2, &c3, &c4);
	      fputc(c3, fout);
	      fputc(c4, fout);
	      len += 2;
	      break;
	    case 0xB0:
	    case 0xC0:
	    case 0xD0:
	    case 0xE0:
	      fprintf(stderr, "ignored non-Ascii character: %2x\n", c1);
	      break;
	    case 0xF0:
	      switch (c1)
		{
		case 0xFD:
		case 0xFE:
		case 0xFF:
		  fprintf(stderr, "ignored non-Ascii character: %2x\n", c1);
		  break;
		default:
		  c2 = fgetc(fin);
		  fprintf(stderr, "ignored user defined SGB code: %2x%2x\n",
			  c1, c2);
		  break;
		}
	    }
	}
#endif
#ifdef DOS
        {
	  if (termStyle)
	    {
	      if (GBmode && len>MAXLEN-5)
		{
		  fprintf(fout, "~}~\n");
		  GBmode = false; len = 0;
		}
	      else if (!GBmode && len>MAXLEN-7)
		{
		  fprintf(fout, "~\n");
		  GBmode = false; len = 0;
		}
	    }
	  if (!GBmode) /* switch to GB mode */
	    {
	      fprintf(fout, "~{");
	      len += 2;
	    }
	  GBmode = true;
	  c2 = fgetc(fin);
	  dos2gb(c1, c2, &c3, &c4);
	  fputc(c3, fout);
	  fputc(c4, fout);
	  len += 2;
	}
#endif
      /* c1 is ASCII */
    else 
	{
	  if (GBmode) {fprintf(fout, "~}"); len += 2;}
	  /* assert(len<=MAXLEN-1) */
	  if (termStyle && (len>MAXLEN-2 || len>MAXLEN-3 && c1=='~'))
	    {
	      fprintf(fout, "~\n");
	      len = 0;
	    }
	  GBmode = false;
	  if (CR2LF && c1=='\r') c1 = '\n';
	  fputc(c1, fout);
	  len++;
	  if (c1=='\n') len=0;
	  else if (c1== '~') {fputc('~', fout); len++;}
	}
    }
  if (GBmode) fprintf(fout, "~}");
}

#ifdef MAC
mac2gb(hi, lo, hi1, lo1)
     int hi, lo, *hi1, *lo1;
{
  if (lo >= 0x9F)
    {
      *hi1 = 0x21 + (hi - 0x81) * 2 + 1;
      *lo1 = 0x21 + (lo - 0x9F);
    }
  else
    {
      *hi1 = 0x21 + (hi - 0x81) * 2;
      if (lo > 0x7F) lo--;
      *lo1 = 0x21 + (lo - 0x40);
    }
}
#endif

#ifdef DOS
dos2gb(hi, lo, hi1, lo1)
     int hi, lo, *hi1, *lo1;
{
  *hi1 = hi - 0x80;
  *lo1 = lo - 0x80;
}
#endif



