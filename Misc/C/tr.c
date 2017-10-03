#include <time.h>
#include <sys/types.h>
#include <stdio.h>

int main()
{
 char buf[100];
 time_t *tloc;
 time_t dd;
 printf("%d\n",time(tloc));
 dd = time(tloc);
 printf("%s",ctime(&dd));
 return(1);
 }
