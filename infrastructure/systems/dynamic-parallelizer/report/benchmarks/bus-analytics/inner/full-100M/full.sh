#!/bin/bash

# full bus analytics script

# 1.sh Vehicles on the road per day
cat $INPUT_FILE |
  sed 's/T..:..:..//' |
  cut -d ',' -f 1,3 |
  sort -u |
  cut -d ',' -f 1 |
  sort |
  uniq -c |
  awk -v OFS="\t" "{print \$2,\$1}"

# 2.sh Days a vehicle is on the road
cat $INPUT_FILE |
  sed 's/T..:..:..//' |
  cut -d ',' -f 3,1 |
  sort -u |
  cut -d ',' -f 2 |
  sort |
  uniq -c |
  sort -k1n |
  awk -v OFS="\t" "{print \$2,\$1}"

# 3.sh Hours each vehicle is on the road
cat $INPUT_FILE |
  sed 's/T\(..\):..:../,\1/' |
  cut -d ',' -f 1,2,4 |
  sort -u |
  cut -d ',' -f 3 |
  sort |
  uniq -c |
  sort -k1n |
  awk -v OFS="\t" "{print \$2,\$1}"

# 4.sh -- Hours monitored each day
cat $INPUT_FILE |
  sed 's/T\(..\):..:../,\1/' |
  cut -d ',' -f 1,2 |
  sort -u |
  cut -d ',' -f 1 |
  sort |
  uniq -c |
  awk -v OFS="\t" "{print \$2,\$1}"

# 5.sh -- Records are day, hour, line, bus
<$INPUT_FILE sed 's/T\(..\):..:../,\1/' | awk -F, '
!seen[$1 $2 $4] { seen[$1 $2 $4] = 1; hours[$1 $4]++; bus[$4] = 1; day[$1] = 1; }
END {
   PROCINPUT_FILEFO["sorted_in"] = "@ind_str_asc"
   for (d in day)
     printf("\t%s", d);
   printf("\n");
   for (b in bus) {
     printf("%s", b);
     for (d in day)
       printf("\t%s", hours[d b]);
     printf("\n");
   }
}'
