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
	ret = chdir(TMPDIR);
	if (ret < 0)
		exit(1);
	fd = syscall(SYS_open, "a", O_RDONLY);
	if (fd < 0)
		exit(1);
	close(fd);
	return 0;
}
