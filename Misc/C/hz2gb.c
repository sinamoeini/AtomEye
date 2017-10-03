static char version[] = "hz2gb 2.0 (July 7, 1992)";

/*
  Copyright (C) 1989, 1992      Fung F. Lee

  hz2gb 2.0: convert a HZ file into a Macintosh* / CCDOS SGB file.
  *For Macintosh pre-6.0.x Simplified Chinese Operating System.
  Later versions use the same internal code (High-bit-set GB) as CCDOS does.

  Usage: hz2gb [-n] [-v] [-e] [-8] [-t]

    -n: convert <LF> to <CR> (from UNIX to Macintosh format)
    -v: verbose mode to report detected errors
    -e: set error reporting channel to stderr instead of stdout
    -8: parse all 8 bits of a character
    -t: ignore line continuation markers (for terminals with limited line size
        and wrap-around capability)

  The HZ specification does not dictate how to convert invalid HZ files,
  just as the definition of a programming language usually does not specify
  how a compiler should handle illegal programming constructs.
  The error recovery procedure of this HZ decoder was designed after
  examination of the conversion errors reported by hz2gb 1.1 of some of the
  "HZ" files posted on the news group alt.chinese.text.  I suspected that 
  most of the errors occured due to improper manual insertion of escape
  sequences, and/or using invalid GB codes, such as those for "space" ($2121).
  Such errors should not have occured if the files were first properly edited
  as GB codes, and then converted by an HZ encoder, such as gb2hz (preferably
  with the -t option.)

  To prevent some hanzi displayers from ill behaviour, the output stream
  should be or should be corrected to be valid mixed ASCII and GB sequences.

  The error recovery procedure is by no means unique, and may change in the
  future. Users should NOT regard the error recovery features as part of the
  HZ specification. 

  This program is free for general distribution.  

  This program runs on UNIX. You are welcome to port it to other operating
  systems.

*/

#include <stdio.h>
#include <string.h>

#define TRUE		1
#define FALSE		0
#define SPACE		0xA1A1		/* GB "space" symbol */
#define BOX		0xA1F5		/* GB "blank box" symbol */
#define isGB1(c)	((c)>=0x21 && (c)<=0x77)	/* GB 1st byte */
#define isGB1U(c)	((c)>=0x78 && (c)<=0x7D)	/* GB 1st byte unused*/
#define isGB2(c)	((c)>=0x21 && (c)<=0x7E)	/* GB 2nd byte */
#define isprint(c)	((c)>=0x20 && (c)<=0x7E)	/* as in <ctype.h> */
#define HI(code)	(((code) & 0xFF00)>>8)
#define LO(code)	((code) & 0x00FF)
#define DB(hi,lo)       (((hi)&0xFF) << 8 | (lo)&0xFF)
#define CLEAN7(c)	((c) & 0x7F)			/* strip MSB */

int LF2CR = FALSE;	/* flag for converting ASCII <LF> to <CR> */
int verbose = FALSE;	/* flag for verbose mode to report errors in input */
int pass8 = FALSE;	/* flat for parsing all 8 bits of a character */
int termStyle = FALSE;	/* flag for ignoring line-continuation markers */
int errorCount = 0;	/* number of parsing errors detected */
FILE *ferr = stdout;	/* error reporting channel */

void usage(), filter(), EOFerror(), ESCerror(), GBerror(), GBerror1();

void usage()
{
    fprintf(stderr, "This is %s. Copyright 1989-1992 Fung F. Lee\n\n", version);
    fprintf(stderr, "usage: hz2gb [-n] [-v] [-e] [-8]\n");
    fprintf(stderr, "-n: convert <LF> to <CR> (from UNIX to Macintosh format)\n");
    fprintf(stderr, "-v: verbose mode to report detected errors\n");
    fprintf(stderr, "-e: set error reporting channel to stderr instead of stdout\n");
    fprintf(stderr, "-8: parse all 8 bits of a character\n");
    fprintf(stderr, "-t: ignore line continuation markers (for terminals \n");
    fprintf(stderr, "    with limited line size and wrap-around capability)\n");
    exit(1); 
}

main(argc, argv)
int argc;
char *argv[];
{
    int i;

    for (i=1; i<argc; i++)
    {
	if (strcmp(argv[i], "-n") == 0)
	    LF2CR = TRUE;
	else if (strcmp(argv[i], "-v") == 0)
	    verbose = TRUE;
	else if (strcmp(argv[i], "-e") == 0)
	    ferr = stderr;
	else if (strcmp(argv[i], "-8") == 0)
	    pass8 = TRUE;
	else if (strcmp(argv[i], "-t") == 0)
	    termStyle = TRUE;
	else
	    usage();
    }
    filter(stdin, stdout);
    if (errorCount > 0)
	fprintf(ferr, "\nNumber of syntax errors caught and recovered = %d\n",
		errorCount);
}

void filter(fin, fout)
FILE *fin, *fout;
{
int c1, c2, c3, c4;
int ASCIImode = TRUE;
    
    while ((c1=fgetc(fin)) != EOF)
    {
	if (!pass8) c1 = CLEAN7(c1);
	if (ASCIImode)
	{
	    if (c1 == '~')
	    {
		if ((c2 = fgetc(fin)) == EOF) {EOFerror(); break;}
		if (!pass8) c2 = CLEAN7(c2);
		switch (c2)
		{
		case '~' : fputc('~', fout); break;
		case '{' : ASCIImode = FALSE; break;
		case '\n': /* line-continuation marker: eat it unless ... */
		    if (termStyle) fputc('\n', fout);
		    break;
		default  : ESCerror(c2);
		    fputc('~', fout); fputc(c2, fout); break;
		}
	    }
	    else
	    {
		if (LF2CR && c1=='\n') c1 = '\r';
		fputc(c1, fout);
	    }
	}
	else /* GBmode */
	{
	    if (isprint(c1))
	    {
		if ((c2 = fgetc(fin)) == EOF) {EOFerror(); break;}
		if (!pass8) c2 = CLEAN7(c2);
		if (isGB1(c1) && isGB2(c2))
		{
		    GBtoSGB(c1, c2, &c3, &c4);
		    fputc(c3, fout);
		    fputc(c4, fout);
		}
		else if (c1 == '~' && c2 == '}')  /* 0x7E7D */
		{
		    ASCIImode = TRUE;
		}  
		else if (isGB1U(c1) && isGB2(c2))  /* 0x78?? - 0x7D?? */
		{
		    GBerror(c1, c2);	/* non-standard extended code? */
		    fputc(HI(BOX), fout); fputc(LO(BOX), fout);
		}
		else if (c1 == '~')	/* 0x7E */
		{
		    GBerror(c1, c2);	/* undefined shift-out code? */
		    ASCIImode = TRUE;	/* safer assumption? */
		    fputc(c1, fout); fputc(c2, fout);
		}
		else if (c1 == ' ')	/* 0x20 */
		{
		    GBerror(c1, c2);	/* looks like artifacts of zwdos? */
		    fputc(c2, fout);
		}
		else if (c2 == ' ')	/* 0x20 */
		{
		    GBerror(c1, c2);	/* null image looks like "sp"? */
		    fputc(HI(SPACE), fout); fputc(LO(SPACE), fout);
		}
		else			/* isprint(c1) && !isprint(c2)) */
		{
		    GBerror(c1, c2);	/* premature shift-out? */
		    ASCIImode = TRUE;	/* safer assumption? */
		    fputc(c1, fout); fputc(c2, fout);
		}
	    }
	    else    /* !isprint(c1) */
	    {
		GBerror1(c1);		/* premature shift-out? */
		ASCIImode = TRUE;	/* safer assumption? */
		fputc(c1, fout);
	    }
	}
    }
}

GBtoSGB(hi, lo, hi1, lo1)
int hi, lo, *hi1, *lo1;
{
#ifdef DOS
    *hi1 = 0x80 | hi;
    *lo1 = 0x80 | lo;
#endif
#ifdef MAC
    *hi1 = 0x81 + (hi - 0x21)/2;
    if (hi%2 != 0)
    {
	*lo1 = 0x40 + (lo - 0x21);
	if (*lo1 >= 0x7F) *lo1 += 1;
    }
    else
	*lo1 = 0x9F + (lo - 0x21);
#endif
}

void EOFerror()
{
    errorCount++;
    if (verbose)
	fprintf(ferr, "\nUnexpected EOF\n");
}

void ESCerror(c)
int c;
{
    errorCount++;
    if (verbose)
	fprintf(ferr, "\nInvalid ASCII escape sequence:\"~%c\"\n", c);
}

void GBerror(c1, c2)
int c1, c2;
{
    errorCount++;
    if (verbose)
	fprintf(ferr, "\nInvalid GB code:\"%c%c\"(0x%4x)\n", c1,c2, DB(c1,c2));
}

void GBerror1(c)
int c;
{
    errorCount++;
    if (verbose)
	fprintf(ferr, "\nInvalid GB code first byte:'%c'(0x%2x)\n", c, c);
}

