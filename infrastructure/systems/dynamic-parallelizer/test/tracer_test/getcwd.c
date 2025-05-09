#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <syscall.h>

#define TMPDIR "/tmp/hs_tracer_test"

void reset(void)
{
	int ret;
	ret = system("rm -rf " TMPDIR);
	ret = system("mkdir " TMPDIR);
}

int main(void)
{
	int fd, ret;
	char name[1024];
	reset();
	ret = chdir(TMPDIR);
	if (ret < 0)
		exit(1);
	ret = syscall(SYS_getcwd, name, 1024);
	if (ret < 0)
		exit(1);
	return 0;
}
