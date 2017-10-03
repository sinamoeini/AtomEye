#include <errno.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

extern int errno;

void main(int argc, char *argv[])
{
    int id ;
    int i;
    struct shmid_ds buf;
    char c;
    int uid, gid;
    int isDelete=0;
    
/* Get args */
    
    while((c=getopt(argc, argv,"hd"))!=(char)(-1))
        switch(c)
        {
        case 'h':
            printf("Version: 1.11, Author: Liao Dongyi(liaody@mit.edu)\n"
                   "  Check and manipulate shared memory blocks\n"
                   "Usage: %s -hd\n"
                   "-h\tThis message\n"
                   "-d\tDelete every found shared memory if possible\n"
                   /*"-c <k>,<s>\tCreate a shared memory of key <k> and size <s>"*/
                   ,argv[0]);
            exit(0);
        case 'd':
            isDelete=1;
            break;
        case 'c':;
        }

/* Initialize, getuid, etc..*/
    printf("User ID=%d, Group ID=%d\n",uid=getuid(), gid=getgid());
/* Try the shmctl(IPC_STAT) to look at all the shared memory */
    printf("****Looking for shared-memory id=0--65535...\n");
    for(i=0;i<65536;i++)
        if(shmctl(i,IPC_STAT,&buf)!=-1)
        {
            printf("  shmctl successed on id=%d.\n", i);
            printf("    Size : %d\n",buf.shm_segsz);
            printf("    Time : %s",ctime(&buf.shm_atime));
            printf("    Key  : %d\n",buf.shm_perm.key);
            printf("    UGID : %d,%d(%s)\n",buf.shm_perm.uid,buf.shm_perm.gid,
                   (uid==buf.shm_perm.uid && gid==buf.shm_perm.gid)?
                   "own":"not own");
            if(isDelete)
                if(shmctl(i,IPC_RMID,NULL)==-1)
                    printf("    - Fail to remove\n");
                else printf("    - Success to remove\n");
        }
        else if(errno != EINVAL && errno != EIDRM)
        {
            printf ("  shmctl failed on id=%d (%s).\n", i,
                    strerror(errno)) ;
            if(isDelete)
                if(shmctl(i,IPC_RMID,NULL)==-1)
                    printf("    - Fail to remove\n");
                else printf("    - Success to remove\n");
        }            
}
