#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define MAXCHAR 50000
 
struct Namelist {
  char * email;
  char * last;
  char * name;
  struct Namelist * next;
  struct Namelist * prev;
};
char buf[MAXCHAR+1], *p, *q;

void insert (Root, Name)
struct Namelist *Root;
struct Namelist *Name;
{
struct Namelist *y = Root->next, *c;

while ((y!=NULL)&&((strcmp(Name->last,y->last)>0)||((strcmp(Name->last,y->last)==0)&&(strcmp(Name->name,y->name)>0))))
{
  c=y;
  y=y->next;
}

if (y!=NULL)
{
   y->prev->next = Name;
   Name -> prev = y->prev;
   Name -> next = y;
   y -> prev = Name;
}
else 
{
  c->next = Name;
  Name -> prev = c;
  Name -> next = NULL;
}

}

main (argc,argv)
int argc;
char **argv;
{
struct Namelist *list=NULL, *a;
char * cc;
FILE *file;

file = fopen(argv[1],"r");
fread (buf,MAXCHAR,1,file);
buf[MAXCHAR] = 0;

list = (struct Namelist *)malloc(sizeof(struct Namelist));
list -> prev = NULL;
list -> next = NULL;
/* the root is an empty node */

q = buf;

while (p=strchr(q,')'))
{
 if (q!=buf) cc = strchr(q,'\n');
  else cc=buf;
 while (!isalpha(*cc)) ++cc;
 
 *p = 0;
 q = p+1;

 while (*(--p)!=' ');
 
 a = (struct Namelist *)malloc(sizeof(struct Namelist));
 a -> last = p+1;
 while (*(--p)!='(');

 *p = 0;
 a -> name = p+1;
 a -> email = cc;

 if (list->next==NULL) 
 {
  list -> next = a;
  list -> next -> prev = list;
  list -> next -> next = NULL;
 }
else insert (list, a);
}

a = list->next;
while (a!=NULL)
{
 printf("\t%-25s %s\n", a->name, a->email);
 a = a->next;
}

}





