window=15

/usr/bin/time -o 100echos2.hs.time -f %e ../../../pash-spec.sh -d2 --window $window 100echos2.sh &> 100echos2.log
/usr/bin/time -o 100echos.hs.time -f %e ../../../pash-spec.sh -d2 --window $window 100echos.sh &> 100echos.log
/usr/bin/time -o giant_file.100.hs.time -f %e ../../../pash-spec.sh -d2 --window $window giant_file.sh &> giant_file.100.log
/usr/bin/time -o giant_file.10000.hs.time -f %e ../../../pash-spec.sh -d2 --window $window giant_file2.sh &> giant_file.1000.log
/usr/bin/time -o multi_files.hs.time -f %e ../../../pash-spec.sh -d2 --window $window 100echos2.sh &> 100echos2.log

/usr/bin/time -o 100echos2.sh.time -f %e sh 100echos2.sh &> 100echos2.log
/usr/bin/time -o 100echos.sh.time -f %e sh 100echos.sh &> 100echos.log
/usr/bin/time -o giant_file.100.sh.time -f %e sh giant_file.sh &> giant_file.100.log
/usr/bin/time -o giant_file.10000.sh.time -f %e sh giant_file2.sh &> giant_file.1000.log
/usr/bin/time -o multi_files.sh.time -f %e sh 100echos2.sh &> 100echos2.log
