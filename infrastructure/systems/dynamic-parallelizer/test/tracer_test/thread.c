#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <syscall.h>
#include <pthread.h>

#define TMPDIR "/tmp/hs_tracer_test"

void reset(void)
{
 	int ret;
	ret = system("rm -rf " TMPDIR);
	ret = system("mkdir " TMPDIR);
	ret = system("mkdir " TMPDIR "/a");
	ret = system("mkdir " TMPDIR "/b");
	ret = system("touch " TMPDIR "/a/f");
}

void *threaded_chdir(void *p)
{
	int ret;
	ret = chdir(TMPDIR "/a");
	return NULL;
}

int main(void)
{
	pthread_t child;
 	int fd, ret;
	reset();
	ret = chdir(TMPDIR "/b");
	ret = pthread_create(&child, NULL, threaded_chdir, NULL);
	ret = pthread_join(child, NULL);
	fd = syscall(SYS_open, "f", O_RDONLY);
	fd = syscall(SYS_open, "g", O_RDONLY);
	if (fd < 0)
                exit(1);
	close(fd);
	return 0;
}
