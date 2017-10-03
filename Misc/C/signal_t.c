/***********************************************/
/* Signal handling in multi-threaded program   */
/* that works for RedHat Linux 5.2 and SGI O2  */
/* % cc signal_t.c -lpthread; a.out            */
/*                                             */
/*                           Ju Li (Mar.2,99)  */
/***********************************************/

#include <stdio.h>
#include <signal.h>
#include <pthread.h>
#define V_SIGNAL SIGUSR1

void nothing_thread (void *a)
{
    printf ("new thread %d created...\n", getpid());
    sleep (1);
    printf ("new thread %d exits.\n\n", getpid());
}

void signal_handler (int i)
{
    pthread_t nothing_thread_ID;
    printf ("\nEntering your own signal_handler\n");
    printf ("you can do some ordinary stuff here,\n");
    printf ("or even create a new thread...\n");
    pthread_create (&nothing_thread_ID, NULL,
		    (void *(*)(void *))nothing_thread,
		    NULL);
}

void signal_handler_thread (void *a)
{
    int sig;
    sigset_t v_signal;
    sigemptyset (&v_signal); 
    sigaddset (&v_signal, V_SIGNAL);
    printf ("Execute \"%% kill -%d %d\" in another shell.\n",
	    V_SIGNAL, getpid());
    while (1)
    {
	sigwait (&v_signal, &sig);
	signal_handler (V_SIGNAL); 
    }
}

void main()
{
    pthread_t signal_handler_thread_ID; 
    printf ("main thread %d comes in\n", getpid());
    pthread_create (&signal_handler_thread_ID, NULL,
		    (void *(*)(void *))signal_handler_thread,
		    NULL);
    printf ("main thread %d comes out\n", getpid());
    while(1);
}
