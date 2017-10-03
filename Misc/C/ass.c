/************************************************/
/* Program Assassin:                            */
/* List and replace hexadecimal bytes of a file */
/*                                              */
/* % cc -O3 -o ass ass.c; ass ass_t 6573 y      */
/*                                              */
/* To check results, use                        */
/* % emacs ass_t --eval "(hexl-mode)" &         */
/*                                              */
/*           Ju Li <liju99@mit.edu> May 4 1999  */ 
/************************************************/

#include <stdio.h>
#include <string.h>

#define ESC	"\x1b"
/*  Foreground Colors  */
#define BLK ESC"[30m"          /* Black    */
#define RED ESC"[31m"          /* Red      */
#define GRN ESC"[32m"          /* Green    */
#define YEL ESC"[33m"          /* Yellow   */
#define BLU ESC"[34m"          /* Blue     */
#define MAG ESC"[35m"          /* Magenta  */
#define CYN ESC"[36m"          /* Cyan     */
#define WHT ESC"[37m"          /* White    */
/*   Hi Intensity Foreground Colors   */
#define HIK ESC"[1;30m"	       /* Black	   */
#define HIR ESC"[1;31m"        /* Red      */
#define HIG ESC"[1;32m"        /* Green    */
#define HIY ESC"[1;33m"        /* Yellow   */
#define HIB ESC"[1;34m"        /* Blue     */
#define HIM ESC"[1;35m"        /* Magenta  */
#define HIC ESC"[1;36m"        /* Cyan     */
#define HIW ESC"[1;37m"        /* White    */
/* High Intensity Background Colors  */
#define HBRED ESC"[41;1m"       /* Red      */
#define HBGRN ESC"[42;1m"       /* Green    */
#define HBYEL ESC"[43;1m"       /* Yellow   */
#define HBBLU ESC"[44;1m"       /* Blue     */
#define HBMAG ESC"[45;1m"       /* Magenta  */
#define HBCYN ESC"[46;1m"       /* Cyan     */
#define HBWHT ESC"[47;1m"       /* White    */
/*  Background Colors  */
#define BBLK ESC"[40m"          /* Black    */
#define BRED ESC"[41m"          /* Red      */
#define BGRN ESC"[42m"          /* Green    */
#define BYEL ESC"[43m"          /* Yellow   */
#define BBLU ESC"[44m"          /* Blue     */
#define BMAG ESC"[45m"          /* Magenta  */
#define BCYN ESC"[46m"          /* Cyan     */
#define BWHT ESC"[47m"          /* White    */
#define NOR ESC"[2;37;0m"      /* Puts everything back to normal */
/* Additional ansi Esc codes added to ansi.h by Gothic  april 23,1993 */
/* Note, these are Esc codes for VT100 terminals, and emmulators */
/* and they may not all work within the mud               */
#define BOLD ESC"[1m"          /* Turn on bold mode */
#define CLR ESC"[2J"           /* Clear the screen */
#define HOME ESC"[H"           /* Send cursor to home position */
#define REF CLR HOME            /* Clear screen and home cursor */
#define BIGTOP ESC"#3"         /* Dbl height characters, top half */
#define BIGBOT ESC"#4"         /* Dbl height characters, bottem half */
#define SAVEC ESC"[s"           /* Save cursor position */
#define REST ESC"[u"            /* Restore cursor to saved position */
#define REVINDEX ESC"M"        /* Scroll screen in opposite direction */
#define SINGW ESC"#5"          /* Normal, single-width characters */
#define DBL ESC"#6"            /* Creates double-width characters */
#define FRTOP ESC"[2;25r"      /* Freeze top line */
#define FRBOT ESC"[1;24r"      /* Freeze bottom line */
#define UNFR ESC"[r"           /* Unfreeze top and bottom lines */
#define BLINK ESC"[5m"         /* Initialize blink mode */
#define REV ESC"[7m"           /* Turns reverse video mode on */
#define HIREV ESC"[1,7m"       /* Hi intensity reverse video */ 
 
#define MAX_PAT 128
#define MAX_STR (MAX_PAT*1024)
/* for debugging */
/* #define MAX_STR (3) */
#define CONTEXT_BYTES 4
#define min(x,y) ((x)-(y)<0?(x):(y))
#define max(x,y) ((x)-(y)>0?(x):(y))

static int hexadecimals_to_bytes
(char *hexadecimals, char *bytes);
static char *bytes_to_hexadecimals
(char *bytes, char *hexadecimals, int n);

/* convert string "3a21" to corresponding 2-byte string */
int hexadecimals_to_bytes (char *hexadecimals, char *bytes)
{
    unsigned int i, n, m;
    char buf[3];
    n = strlen (hexadecimals);
    if ( n & 1 )
    {
	printf ("Hexadecimal string \"%s\" should have even length, "
		"e.g. 1F2A.\n", hexadecimals); 
	exit(1);
    }
    n /= 2;
    for (i=0; i<n; i++)
    {
        buf[0] = hexadecimals[2*i];
        buf[1] = hexadecimals[2*i+1];
        buf[2] = '\0';
	sscanf (buf, "%x", &m);
	bytes[i] = m;
    }
    bytes[n] = '\0';
    return (n);
} /* end hexadecimals_to_bytes() */


/* convert 2-byte string to "3a21" */
char *bytes_to_hexadecimals (char *bytes, char *hexadecimals, int n)
{
    int i;
    for (i=0; i<n; i++)
	sprintf (hexadecimals+2*i, "%02x", ((int)bytes[i])&0xff);
    hexadecimals[2*n] = '\0';
    return(hexadecimals);
} /* end bytes_to_hexadecimals() */


#ifndef NO_ASS_TEST
int main (int argc, char *argv[])
{
    int i, j, n, m, coincidence;
    char pat[MAX_PAT], buf[MAX_STR], fn_out[128], tmp[2*MAX_PAT],
	tmp_pat[MAX_PAT], *p, *q, *v, *w, *f;
    FILE *in, *out;
    
    if (argc < 3)
    {
        printf ("Assassinate a sequence of hexadecimal bytes "
                "in a binary file.\n");
        printf ("Usage: %% %s file_in CAC7D2E4 [file_out]\n", argv[0]);
	exit(1);
    }
    else if (argc == 3) /* default output filename */
	sprintf (fn_out, "%s.a", argv[1]);
    else if ( !strcmp(argv[1], argv[3]) )
    {  /* cannot read and write at the same time */
	printf ("Input and output files must have different names.\n",
		argv[1], argv[3]);
	exit(1);
    }
    else sprintf (fn_out, "%s", argv[3]);
    
    if ( (in=fopen(argv[1],"r")) == NULL )
    {
	printf ("File \"%s\" does not exist.\n", argv[1]);
	exit(1);
    }
    out = fopen(fn_out, "w");
    
    /* convert input pattern to bytes format */
    n = hexadecimals_to_bytes (argv[2], pat);
    if (n >= MAX_PAT)
    {
	printf ("Hexadecimal string >= %d bytes, too long.\n",
		MAX_PAT);
	exit(1);
    }
    
    for (coincidence=0, p=buf;;)
    {  /* fill up the buffer with data */
	m = fread (p, sizeof(char), buf+MAX_STR-p, in);
	f = p + m; /* final pointer position */
	for (v=buf; v<f-n; v++)
	    if (memcmp(v, pat, n) == 0)
	    {
		++coincidence;
		printf ("\nFound pattern at byte position %d ("
			HIW"0x%x"NOR"):\n",
			ftell(in)-(f-v), (unsigned)ftell(in)-(f-v));
		w = max(v-CONTEXT_BYTES, buf);
		bytes_to_hexadecimals(w, tmp, v-w);
		printf ("%s  ", bytes_to_hexadecimals(w, tmp, v-w));
		printf (HIG"%s"NOR "("HIC, bytes_to_hexadecimals(v, tmp, n));
		for (i=0; i<n; i++)
		    putchar(isprint(pat[i])?pat[i]:'.');
		printf (NOR ")  ");
		w = min(v+n+CONTEXT_BYTES, f);
		printf ("%s\n", bytes_to_hexadecimals(v+n, tmp, w-v-n));
		while (1)
		{
		    printf ("Replace by (press return to skip):\n");
		    for (i=0; (tmp[i]=getc(stdin))!='\n'; i++); tmp[i]='\0';
		    /* return pressed */
		    if (tmp[0] == '\0') break;
		    j = hexadecimals_to_bytes (tmp, tmp_pat);
		    if ( j != n )
			printf ("You are not allowed to change the length of "
				"byte-string %s.\n", argv[2]);
		    else
		    {
			memcpy (v, tmp_pat, n);
			v += n-1;
			break;
		    }
		}
	    }
	if ( !feof(in) )
	{ 
	    fwrite (buf, sizeof(char), v-buf, out);
	    fflush(out);
	    /* move the rest of buffer to head */
	    for (p=buf,q=v; q<f; p++,q++) *p=*q;
	}
	else
	{
	    fwrite (buf, sizeof(char), f-buf, out);
	    fclose (in); fclose (out);
	    printf ("\nSaved to file \"%s\".\n", fn_out);
	    break;
	}    
    }
    return (0);
} /* end main() */
#endif
