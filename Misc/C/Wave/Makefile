#C++ Project Makefile  Jan.23,1995, LiJu.

CC      = gcc
LIBS    = -L/mit/recipes/sun4lib -lg++ -lm 
OBJS    = complex.o optical.o

run:		$(OBJS)
	$(CC) -o $@ $(OBJS) $(LIBS) -O4

complex.o:	complex.cc complex.h
	$(CC) -c complex.cc

optical.o:	optical.cc complex.h
	$(CC) -c optical.cc

