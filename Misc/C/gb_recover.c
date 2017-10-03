/*******************************************/
/* gb_recover.c:                           */
/*                                         */
/* Empirical attempt to recover 8-bit mail */
/* delivered by standard 7-bit service.    */
/*******************************************/

#include <stdio.h>

int p;
char stack[2];
int just_gb;


/* The two-byte sequence is a gb encoding or not */
int common_gb (char *stack)
{
    unsigned int i;
    i = ((stack[0]|0x80)<<8) + (stack[1]|0x80);
    return ( (i==0xA1A3) ||   /* . */ 
             (i==0xA1B0) ||   /* " */
             (i==0xA1B1) ||   /* " */
             (i==0xA3A1) ||   /* ! */ 
             (i==0xA3A8) ||   /* ( */
             (i==0xA3A9) ||   /* ) */
             (i==0xA3AC) ||   /* , */
             (i==0xA3BA) ||   /* : */
             (i==0xA3BB) ||   /* ; */
             (i==0xA3BF) ||   /* ? */
             ((i>0xB0A1) && (i<0xD7FE)) );
} /* end common_gb() */


/* When the 2-byte stack is full... */
void handle_fullstack()
{
    if ( (stack[0]=='\\') && (stack[1]=='n') )
    {
        putc('\n',stdout);
        p = 0;
        just_gb = 0;
    }
    else if (common_gb(stack))
    {
        putc(stack[0]|0x80,stdout);
        putc(stack[1]|0x80,stdout);
        p = 0;
        just_gb = 1;
    }
    else
    {
        putc(stack[0],stdout);
        stack[0] = stack[1];
        p--;
        just_gb = 0;
    }
    return;
} /* end handle_fullstack() */


int main (int argc, char *argv[])
{
    FILE *fp;
    int i,j;
    if (argc==2) fp = fopen(argv[1],"r");
    else fp = stdin;
    p = 0;
    just_gb = 0;
    while ((i=fgetc(fp)) != EOF)
    {
        if ( (!just_gb) && (i=='>') )
        {
            if ( (j=fgetc(fp)) == '^')
            {
                if (p==1)
                {
                    putc(stack[0],stdout);
                    p = 0;
                }
                putc(i|0x80,stdout);
                putc(j|0x80,stdout);
                just_gb = 1;
            }
            else
            {
                ungetc(j,fp);
                if (p==1)
                {
                    putc(stack[0],stdout);
                    p = 0;
                }
                putc(i,stdout);
                just_gb = 0;
            }
        }
        else
        {
            stack[p++] = i;
            if (p==2) handle_fullstack();
        }
    }
    return (0);
} /* end main() */

/* cc gb_recover.c -o ~/Co/linuxBin/gb_recover */
