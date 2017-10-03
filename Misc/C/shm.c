#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>
#include <string.h>
#define KEY 329

/* define common variables */
unsigned int ID;
double *shm_block;
double *u,*u2,*u4,*u6,*u8,*u10,*uold,*copy;
int size_of_vector_in_bytes,size_of_vector_in_double;
int *schedule_flags;
#define DO_NOTHING       0
#define U_TO_MY_COPY     1
#define U2_TO_MY_COPY    2
#define U4_TO_MY_COPY    3
#define U6_TO_MY_COPY    4
#define U8_TO_MY_COPY    5
#define SLAVE_DIE       99

extern struct
{
    unsigned int SU,SU2,SU4,SU6,SU8,SU10,SUOLD,SCOPY;
} offset_ ;

extern struct
{
    int num_processors;
    int *domain;
} layout_ ;

void Error(const char *fmt,...)
{
    printf("Error:%s(%s)\n",fmt,strerror(errno));
}

int generate_shm_block_ (int *n3, int *num_processors)
{
    unsigned int size;
    size = sizeof(double)*2*(*n3)*(7+(*num_processors))
	   +sizeof(int)*(*num_processors);
    ID=shmget(KEY,size,IPC_CREAT|0600);
    if(ID!=-1)
    {
	/* system will find a good logical address 
	   to attach the shared memory segment */
	shm_block = shmat(ID,(char *)0,0600);
        if(((int)shm_block)==-1)
	{
            Error("Failed to attach shared memory");
	    return (0);
	}
	else
	{
	    uold = shm_block; 
	    u = uold + 2*(*n3);
	    u2  = u  + 2*(*n3);
	    u4  = u2 + 2*(*n3); 
	    u6  = u4 + 2*(*n3); 
	    u8  = u6 + 2*(*n3); 
	    u10 = u8 + 2*(*n3);
	    copy = u10 + 2*(*n3);
	    schedule_flags = (int *)(copy+2*(*n3)*(*num_processors));
            /* offset list */
	    offset_.SUOLD = 0;
	    offset_.SU = offset_.SUOLD + (*n3);
	    offset_.SU2  = offset_.SU  + (*n3);
	    offset_.SU4  = offset_.SU2 + (*n3);
	    offset_.SU6  = offset_.SU4 + (*n3);
	    offset_.SU8  = offset_.SU6 + (*n3);
	    offset_.SU10 = offset_.SU8 + (*n3);
	    offset_.COPY = offset_.SU10 + (*n3);
	    size_of_vector_in_double = 2*(*n3);
	    size_of_vector_in_bytes  = 2*(*n3)*sizeof(double);
	    return (1);
	}
    }
    else
        Error("Failed to get ** bytes of shared memory");
}

void slave_wait_for_command_ (int *My_idx)
{
    double *in, *out;
    unsigned int start, end;
    while (1)
    {
	switch(schedule_flags[*My_idx-1])
	{
	case U_TO_MY_COPY:   in=u;  goto 21;
	case U2_TO_MY_COPY:  in=u2; goto 21;
	case U4_TO_MY_COPY:  in=u4; goto 21;
	case U6_TO_MY_COPY:  in=u6; goto 21;
        case U8_TO_MY_COPY:  in=u8; 
	21: out = copy+(*My_idx-1)*size_of_vector_in_double;
	f77_sparse_multiply_ (in,out,&layout_.domain[*My_idx-1],
			             &layout_.domain[*My_idx]);
	/* negate schedule_flag to show that it's done. */
	schedule_flags[*My_idx-1] = -schedule_flags[*My_idx-1];
	break;
        case SLAVE_DIE: return;
	}
    }
}
    
void quick_sum_copy (double *add_to_where)
{
    register double a;
    register int jump = size_of_vector_in_double;
    register int i,j,nc = layout_.num_processors;

    for (i=0;i<jump;i++)
    {
	a = copy[i];
	for (j=1;j<nc;j++) a+=copy[i+j*jump];
	add_to_where[i] = a;
    }
}

void master_multiply_all_
(double *v,double *v2,double *v4,double *v6,double *v8,double *v10, int *ibb,
 double *v2_ibb,double *v4_ibb,double *v6_ibb,double *v8_ibb,double *v10_ibb)
{
    /* first copy u to shared memory */
    memcpy(u,   v,   size_of_vector_in_bytes);
    /* zero the higher order u's */
    memset(u2,  0, 5*size_of_vector_in_bytes);
    /* add in perturbation */
    call write_shm_ (v2_ibb, SU2+ibb);
    call write_shm_ (v4_ibb, SU4+ibb);
    call write_shm_ (v6_ibb, SU6+ibb);
    call write_shm_ (v8_ibb, SU8+ibb);
    call write_shm_ (v10_ibb, SU10+ibb);
    /* do the multiplication on many processors */
    call master_multiply(U_TO_MY_COPY);
    call master_multiply(U2_TO_MY_COPY);
    call master_multiply(U4_TO_MY_COPY);
    call master_multiply(U6_TO_MY_COPY);
    call master_multiply(U8_TO_MY_COPY);
    /* copy the results back to single processor */
    memcpy(v2,   u2,   size_of_vector_in_bytes);   
    memcpy(v4,   u4,   size_of_vector_in_bytes);   
    memcpy(v6,   u6,   size_of_vector_in_bytes);   
    memcpy(v8,   u8,   size_of_vector_in_bytes);   
    memcpy(v10,  u10,  size_of_vector_in_bytes);   
}    
     
void master_multiply (int what)
{
    int i,sum;
    /* send multiply command to all slaves */
    mmset(schedule_flags,what,layout_.num_processors);
    /* does his own share (good master) */
    switch(what)
    {
    case U_TO_MY_COPY:   in=u;  goto 21;
    case U2_TO_MY_COPY:  in=u2; goto 21;
    case U4_TO_MY_COPY:  in=u4; goto 21;
    case U6_TO_MY_COPY:  in=u6; goto 21;
    case U8_TO_MY_COPY:  in=u8; 
    21:   f77_sparse_multiply_ (in,copy,layout_.domain,
			                layout_.domain+1);
    break;
    }
    /* check how the slaves are doing */
    while(1)
    {
	sum=0;
	for (i=1;i<layout_.num_processors)
	    {
		sum+=schedule_flags[i];
	    }
	if (sum==(-what)*(layout_.num_processors-1)) break;
    }
    /* sum up the copies */
    switch(what)
    {
    case U_TO_MY_COPY:   quick_sum_copy(u2); break;
    case U2_TO_MY_COPY:  quick_sum_copy(u4); break;
    case U4_TO_MY_COPY:  quick_sum_copy(u6); break;
    case U6_TO_MY_COPY:  quick_sum_copy(u8); break;
    case U8_TO_MY_COPY:  quick_sum_copy(u10); break;
    }
    return;
}

void read_shm_ (double *complex, unsigned int *offset)
{ /* read in a complex number from shared memory */
    *complex = shm_block[(*offset-1)*2];
    *(complex+1) = shm_block[(*offset-1)*2+1];
}

void write_shm_ (double *complex, unsigned int *offset)
{ /* write a complex number into shared memory */
    shm_block[(*offset-1)*2] = *complex;
    shm_block[(*offset-1)*2+1] = *(complex+1);
}

void zero_shm_block_ (unsigned int *offset, unsigned int *size)
{ /* write a complex number into shared memory */
    memset(shm_block+(*offset)*2,0,(*size)*sizeof(double)*2);
}

void kill_all_slaves_ ()
{
    mmset(schedule_flags,SLAVE_DIE,layout_.num_processors);
}

void free_shm_block_()
{
    shmctl (ID,IPC_RMID,NULL);
    return;
}

