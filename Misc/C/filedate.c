/*********************************************************************/
/* Print the last modification date of a file in format "2006-12-14" */
/*                                                                   */
/* cc -static -o filedate filedate.c; filedate filedate.c            */
/*                                                                   */
/*                             Ju Li <li.562@osu.edu> Dec. 14, 2006  */ 
/*********************************************************************/

#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>

int main (int argc, char *argv[])
{
    char *date; 
    int res; 
    struct stat fstatdata; 
    int fh; 
    struct tm* tmdata;
    
    if (argc != 2)
    {
        printf ("Print file last modification date in format "
                "\"2006-12-14\".\n");
        printf ("Usage: %% %s file_in\n", argv[0]);
	exit(1);
    }

    if ((fh = open (argv[1], O_RDONLY)) < 0) 
    {   
        perror (" Error opening file"); 
        exit (EXIT_FAILURE); 
    } 
    
    if ((res = fstat (fh, &fstatdata)) != 0) 
    {   
        perror (" Failure calling fstat"); 
        exit (EXIT_FAILURE); 
    } 

    tmdata = localtime (&fstatdata.st_mtime);
    printf ("%04d-%02d-%02d\n",
            tmdata->tm_year+1900, tmdata->tm_mon+1, tmdata->tm_mday);
    
    /* date = asctime (localtime (& fstatdata.st_ctime));  */
    /* printf ("\nDate: %s", date);  */
    /* printf (" Mode: %d\n", fstatdata.st_mode);  */
    /* printf (" Size: %ld\n", fstatdata.st_size);  */
    
    return (0);
} /* end main() */
