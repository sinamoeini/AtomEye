/*
** A program to facilitate file operations.
*/

#include <stdio.h>
#include <stdlib.h>

/*
** fileread() reads the file with Filename into the Buffer.
** Buffer is automatically allocated memory.(REMEMBER to free!)
** Return 1 if success.
*/
int fileread(char * Filename, char * Buffer);
int filewrite(char * Buffer, char * Filename);

int fileread(char * Filename, char * Buffer)
{
  FILE *file;
  long size;
  file = fopen(Filename,"r");
  if (file==NULL)
    {
      fprintf(stderr,
	      "Error: cannot locate file %s\n.",
	      file);
      return(0);
    }
  fseek(file,0L,SEEK_END);
  size = ftell(file);
  Buffer = (char *)malloc(size);
  if (Buffer==NULL) 
    {
      fprintf(stderr,
	      "Error: don't have enough memory to load in %s\n.",
	      file);
      fclose(file);
      return(0);
    }
  rewind(file);
  fread(Buffer,size,1,file);
  fclose(file);
  return(1);
}


int filewrite(char * Buffer, char * Filename)
{
  FILE *file;
  long l,size=0;
  file = fopen(Filename,"w");
  while (Buffer[size]!=EOF) size++;
  l=fwrite(Buffer,1,size,file);
  if (l!=size) 
    {
      fprintf(stderr,
	      "Error: failed to write in all the file %s\n.",
	      file);
      fclose(file);
      return(0);
    }
  fclose(file);
  return(1);
}




