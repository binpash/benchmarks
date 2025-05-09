#include <stdio.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/time.h>

#define NUM_CALLS 100000

int main() {
    struct timeval start, end;
    long seconds, useconds;
    double mtime;

    gettimeofday(&start, NULL); // get the start time
    for (int i = 0; i < NUM_CALLS; ++i) {
        syscall(SYS_getpid);
    }
    gettimeofday(&end, NULL);  // get the end time

    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

    printf("Elapsed time for getpid syscall: %.3f milliseconds\n", mtime);
    printf("Average time per getpid syscall: %.3f microseconds\n", mtime * 1000 / NUM_CALLS);
    return 0;
}
