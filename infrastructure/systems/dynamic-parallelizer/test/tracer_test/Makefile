objs = chdir pid fork symlink getcwd retval mkdir renameat2 symlink_syscall

CFLAGS=-O2

.PHONY: all clean

all: $(objs) thread

thread: thread.c
	$(CC) $(CFLAGS) -pthread $< -o $@

clean:
	rm $(objs)
