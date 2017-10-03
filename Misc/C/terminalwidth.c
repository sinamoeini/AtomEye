/*********************************************************/
/* Print the terminal width; if no terminal, status = 1  */
/*                                                       */
/* cc -static terminalwidth.c -o $BIN_PATH/terminalwidth */
/*********************************************************/

#include <stdio.h>
#include <err.h>
#include <sys/ioctl.h>

int main(void)
{
    struct winsize ws;

    if (ioctl(1, TIOCGWINSZ, (void *)&ws) == 0) 
        printf("%d\n", ws.ws_col);
    else
        err(1, "ioctl failed");
    return 0;
}
