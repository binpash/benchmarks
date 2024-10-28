#!/bin/sh

# Move to the src directory
cd /root/redis/src

# Regenerate the release.h file
GIT_SHA1=`(git show-ref --head --hash=8 2> /dev/null || echo 00000000) | head -n1`
GIT_DIRTY=`git diff --no-ext-diff 2> /dev/null | wc -l`
BUILD_ID=`uname -n`"-"`date +%s`
echo "#define REDIS_GIT_SHA1 \"$GIT_SHA1\"" > release.h
echo "#define REDIS_GIT_DIRTY \"$GIT_DIRTY\"" >> release.h
echo "#define REDIS_BUILD_ID \"$BUILD_ID\"" >> release.h

CFLAGS="-g -ggdb -pedantic -DREDIS_STATIC= -std=c11 -Wall -W -Wno-missing-field-initializers -O2 -g -ggdb -I../deps/hiredis -I../deps/linenoise -I../deps/lua/src -I../deps/hdr_histogram -DUSE_JEMALLOC -I../deps/jemalloc/include"
LDFLAGS="-rdynamic ../deps/hiredis/libhiredis.a ../deps/lua/src/liblua.a ../deps/hdr_histogram/hdr_histogram.o ../deps/linenoise/linenoise.o ../deps/jemalloc/lib/libjemalloc.a -lm -lpthread -ldl"

# Build all source files except the ae_* files
SRC="acl.c adlist.c ae.c anet.c aof.c bio.c bitops.c blocked.c childinfo.c cli_common.c cluster.c config.c connection.c crc16.c crc64.c crcspeed.c db.c debug.c defrag.c dict.c endianconv.c evict.c expire.c geo.c geohash.c geohash_helper.c gopher.c hyperloglog.c intset.c latency.c lazyfree.c listpack.c localtime.c lolwut.c lolwut5.c lolwut6.c lzf_c.c lzf_d.c memtest.c module.c monotonic.c mt19937-64.c multi.c networking.c notify.c object.c pqsort.c pubsub.c quicklist.c rand.c rax.c rdb.c redis-benchmark.c redis-check-aof.c redis-check-rdb.c redis-cli.c release.c replication.c rio.c scripting.c sds.c sentinel.c server.c setcpuaffinity.c setproctitle.c sha1.c sha256.c siphash.c slowlog.c sort.c sparkline.c syncio.c t_hash.c t_list.c t_set.c t_stream.c t_string.c t_zset.c timeout.c tls.c tracking.c util.c ziplist.c zipmap.c zmalloc.c"
gcc $CFLAGS -c $SRC

# Build redis-server
OBJ="adlist.o quicklist.o ae.o anet.o dict.o server.o sds.o zmalloc.o lzf_c.o lzf_d.o pqsort.o zipmap.o sha1.o ziplist.o release.o networking.o util.o object.o db.o replication.o rdb.o t_string.o t_list.o t_set.o t_zset.o t_hash.o config.o aof.o pubsub.o multi.o debug.o sort.o intset.o syncio.o cluster.o crc16.o endianconv.o slowlog.o scripting.o bio.o rio.o rand.o memtest.o crcspeed.o crc64.o bitops.o sentinel.o notify.o setproctitle.o blocked.o hyperloglog.o latency.o sparkline.o redis-check-rdb.o redis-check-aof.o geo.o lazyfree.o module.o evict.o expire.o geohash.o geohash_helper.o childinfo.o defrag.o siphash.o rax.o t_stream.o listpack.o localtime.o lolwut.o lolwut5.o lolwut6.o acl.o gopher.o tracking.o connection.o tls.o sha256.o timeout.o setcpuaffinity.o monotonic.o mt19937-64.o"
gcc $CFLAGS -o redis-server $OBJ $LDFLAGS

# Build redis-cli
OBJ="anet.o adlist.o dict.o redis-cli.o zmalloc.o release.o ae.o crcspeed.o crc64.o siphash.o crc16.o monotonic.o cli_common.o mt19937-64.o"
gcc $CFLAGS -o redis-cli $OBJ $LDFLAGS

# Build redis-benchmark
OBJ="ae.o anet.o redis-benchmark.o adlist.o dict.o zmalloc.o release.o crcspeed.o crc64.o siphash.o crc16.o monotonic.o cli_common.o mt19937-64.o"
gcc $CFLAGS -o redis-benchmark $OBJ $LDFLAGS

install redis-server redis-sentinel
install redis-server redis-check-aof
install redis-server redis-check-rdb
