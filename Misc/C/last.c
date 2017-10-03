#include <stdio.h>
#define MAXCHAR 1000000

main(int argc, char *argv[])
{
  FILE * file;
  char buf[MAXCHAR], *p;
  long  size;

  if (argv[1]!=NULL) file=fopen(argv[1],"r");
    else {
      fprintf(stderr, "Must specify filename.\n");
      exit(1);
    }

  size = fread(buf,1,MAXCHAR,file);
  fclose (file);

  if (size==MAXCHAR)
    {
      fprintf(stderr, "Unable to handle such a big file.\n");
      exit(1);
    }

  p=buf+size;

  while (p>buf)
    {
      if (*p=='B')
	{
	  p+=2;
	  break;
	}
      p--;
    }

  file=fopen(argv[1],"w");
  fwrite (p,buf+size-p,1,file);
  fclose (file);
 
}
  
  
  
