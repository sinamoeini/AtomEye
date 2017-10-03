#include <errno.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#define  NULL        ((void *) 0)

typedef  enum {
    True  = 1,
    False = 0
} bool ;

#define  IPCKEY      2347
#define  IPCSIZE     sizeof(INFO)

typedef  int  IPCID ;

typedef  struct {
    int  socket ;
    char whatever [64] ;
} INFO ;

extern   int  errno ;

main ()
{
    IPCID  id ;
    INFO * infoptr ;
    bool   owner  = False ;     /* Responsible to creat and destroy the mem. */
    bool   reader = True ;      /* True for reading, False for writing. */
    
    id = shmget (IPCKEY, IPCSIZE, IPC_CREAT | IPC_EXCL | 0600) ;
    if (id == -1)
        if (errno != EEXIST)
            {
            printf ("shmget failed to create a shared mem, errno = %d.\n", errno) ;
            exit (-1) ;
            }
        else
            {
            id = shmget (IPCKEY, IPCSIZE, 0600) ;
            if (id == -1)
                {
                printf ("shmget failed to get an id, errno = %d.\n", errno) ;
                exit (-3) ;
                }
            }
    else
        {
        owner  = True ;
        reader = False ;
        }
    
    
    infoptr = (INFO *) shmat (id, NULL, SHM_RND) ;
    if (infoptr == NULL || infoptr == (void *) -1)
        {
        printf ("shmat failed, errno = %d.\n", errno) ;
        exit (-3) ;
        }    

    if (reader == False)
        {
        int dummy ;
        infoptr->socket = 9999 ;
        printf ("Please give a dumb number to end the program and move the mem : ") ;
        scanf ("%d", &dummy) ;
        infoptr->socket = dummy ;
        }
    else
        {
       /*
        * This is for demonstration only. The real projects involve 
        * semaphore or at least pipe, which complicates this little 
        * demo so I omit it. You may find info about semaphore in any 
        * Unix textbook, or you can send me e-mail and I'll give you
        * my CSync class for synchronization.
        *
        */
    
        printf ("The socket number is %d.\n", infoptr->socket) ;
        } 
    
    shmdt (infoptr); 

   /*
    * One difficulty in real application of shared memory is who and when
    * to remove the shared memory. Here I simplify the problem by using
    * "owner" variable and persuming that the owner process is executed
    * first and terminated last, which is, in most cases, not true.  So 
    * maybe make a process solely to create and destroy the shared mem
    * is a good idea. Other process just shmget the id without creating
    * it. Moreover, I think there are some other better approaches, let
    * your brain run ...
    *
    * Another issue I do not handle well in this demo program is signal.
    * You have to catch every possible signal to remove the memory, or
    * you may fail to create the memory with same key next time.
    *
    */

    if (owner == True)
        shmctl (id, IPC_RMID) ;
    exit (0) ;
}







