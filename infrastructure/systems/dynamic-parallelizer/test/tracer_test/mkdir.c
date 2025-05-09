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
	reset();
	ret = system("mkdir -p " TMPDIR "/a");
	ret = system("mkdir -p " TMPDIR "/a");
	return 0;
}
