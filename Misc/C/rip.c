/***********************************************/
/* get rid of lines that begins with a numeric */
/***********************************************/
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>

int main(int argc, char *argv[])
{
    char buffer[300];
    int i, math;
    time_t tloc, tt;
    char *dat="dat.out", *rep="rep";
    char *head=" \n \n \n     Simulation Report of Ju Li on ";
    FILE *data, *report;
    if (argc>1)  rep = argv[1];
    if (argc==3) dat = argv[2];
    data = fopen(dat,"w");
    report = fopen(rep,"w");
    fwrite (head,strlen(head),1,report);
    tt = time(&tloc);
    fwrite (ctime(&tt),26,1,report);
    i=0;
    while(!feof(stdin))
    {
	while((buffer[i]=getc(stdin))==' ') i++;
	math=(buffer[i]<='9')&&(buffer[i]>='0');
	while(!(feof(stdin)||(buffer[i]=='\n'))) buffer[++i]=getc(stdin);
	if (feof(stdin)) i--;
	if (math) fwrite(buffer,i+1,1,data);
	else fwrite(buffer,i+1,1,report);
	i=0;
    }
    fclose(data);
    fclose(report);
    return(0);
}   

 
 




