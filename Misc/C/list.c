/*
** list.c 
**
** Li Ju. Jul.10, 95.
** Using sift method to give a prime number list.
*/

#include <stdio.h>
#define MAXPRIME 1000000

#define getbit(p) ( land[p>>3] & (0x1 << (p&7)) )
#define zerobit(p) { land[p>>3] &= 255-(0x1 << (p&7)); } 

unsigned short land[MAXPRIME/8+4];

void main()
{
  long int i,j,k;
  int number[12]={0};

  k = MAXPRIME/8+4;
  for (i=0;i<k;i++)
    land[i] = 255;

  for (i=2;i<MAXPRIME;i++)
    {
      if (getbit(i))
	{
	  for (j=2*i;j<MAXPRIME;j+=i) 
	    zerobit(j);
	  number[i%12]++;
	}
    }
  
  printf ("Small = %d\n", number[5]+number[7]);
  printf ("Large = %d\n", number[1]+number[11]);
  
}



