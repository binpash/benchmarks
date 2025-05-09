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
	ret = system("mkdir " TMPDIR "/b");
	ret = system("touch " TMPDIR "/b/f");
}

int main(void)
{
 	int fd, ret;
	reset();
	/* if (ret == 0) { */
	/* 	ret = chdir(TMPDIR "/a"); */
	/* 	if (ret < 0) */
	/* 		exit(1); */
	/* } else { */
	/* 	ret = chdir(TMPDIR "/b"); */
	/* 	if (ret < 0) */
	/* 		exit(1); */
	/* } */
	ret = chdir(TMPDIR "/a");
	fd = syscall(SYS_symlink, TMPDIR "/b/f", TMPDIR "/a/f");
	if (fd < 0)
                exit(1);
	close(fd);
	return 0;
}
