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
	ret = system("mkdir " TMPDIR "/a");
	ret = system("touch " TMPDIR "/a/f");
}

int main(void)
{
 	int fd, ret;
	reset();
	syscall(SYS_renameat2,
		AT_FDCWD, TMPDIR "/a/f",
		AT_FDCWD, TMPDIR "/a/d", 0);
	return 0;
}
