#!/bin/bash
IN=$IN
IN1=$IN/1.txt
IN2=$IN/2.txt
IN3=$IN/3.txt
IN4=$IN/4.txt
IN5=$IN/5.txt
IN6=$IN/6.txt
IN7=$IN/7.txt
IN8=$IN/8.txt
IN91=$IN/9.1.txt
IN92=$IN/9.2.txt
IN93=$IN/9.3.txt
IN94=$IN/9.4.txt
IN95=$IN/9.5.txt
IN96=$IN/9.6.txt
IN97=$IN/9.7.txt
IN98=$IN/9.8.txt
IN99=$IN/9.9.txt
IN10=$IN/10.txt
IN11=$IN/11.txt
IN12=$IN/12.txt

#1.sh 1.0: extract the last name
cat $IN1 | cut -d ' ' -f 2

#2.sh 1.1: extract names and sort
cat $IN1 | cut -d ' ' -f 2 | sort

#3.sh 1.2: extract names and sort
cat $IN1 | head -n 2 | cut -d ' ' -f 2

#4.sh 1.3: sort top first names
cat $IN1 | cut -d ' ' -f 1 | sort | uniq -c | sort -r

#5.sh 2.1: get all Unix utilities
cat $IN2 | cut -d ' ' -f 4 | tr -d ','

#6.sh 3.1: get lowercase first letter of last names (awk)
cat $IN3 | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

#7.sh 4.1: find number of rounds
cat $IN4 | tr ' ' '\n' | grep '\.' | wc -l

#8.sh 4.2: find pieces captured by Belle
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l

#9.sh 4.3: find pieces captured by Belle with a pawn
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]' | wc -l

#10.sh 4.4: histogram of Belle's captures (-pawns) by each type of piece
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep '[KQRBN]' | cut -c 1-1 | sort | uniq -c | sort -nr

#11.sh 4.5: 4.4 + pawns
cat $IN4 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq -c | sort -nr

#12.sh 4.6: piece used the most by Belle
cat $IN4 | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort -r | uniq | head -n 3 | tail -n 1

#13.sh 5.1: extract hello world
cat $IN5 | grep 'print' | cut -d "\"" -f 2 | cut -c 1-12

#14.sh 6.1: order the bodies by how easy it would be to land on them in Thompson's Space Travel game when playing at the highest simulation scale
cat $IN6 | awk "{print \$2, \$0}" | sort -nr | cut -d ' ' -f 2

#15.sh 7.1: identify number of AT&T unix versions
cat $IN7 | cut -f 1 | grep 'AT&T' | wc -l

#16.sh 7.2: find  most frequently occurring machine
cat $IN7 | cut -f 2 | sort -n | uniq -c | sort -nr | head -n 1 | tr -s ' ' '\n' | tail -n 1

#17.sh 7.3: all the decades in which a unix version was released
cat $IN7 | cut -f 4 | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/

#18.sh 8.1: count unix birth-year
cat $IN8 | tr ' ' '\n' | grep 1969 | wc -l

#19.sh 8.2: find Bell Labs location where Dennis Ritchie had his office
cat $IN8 | grep 'Bell' | awk 'length <= 45' | cut -d ',' -f 2 | awk "{\$1=\$1};1"

#20.sh 8.3: find names of the four people most involved with unix
cat $IN8 | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1 | head -n 1

#21.sh 8.4: find longest words without hyphens
cat $IN8 | tr -c "[a-z][A-Z]" '\n' | sort | awk "length >= 16"

#22.sh # 8.5: Find second-most-freq 8-character word(s) without hyphens
cat $IN8 | tr -cs '[:alpha:]' '\n' | awk 'length($0) == 8' | tr '[:upper:]' '[:lower:]' | sort | uniq -c | sort -nr | awk 'NR==2 {print $2}'

#23.sh 9.1: extract the word PORT
cat $IN91 | tr ' ' '\n' | grep '[A-Z]' | tr '[a-z]' '\n' | grep '[A-Z]' | tr -d '\n' | cut -c 1-4

#24.sh 9.2: extract the word BELL
cat $IN92 | cut -c 1-1 | tr -d '\n'

#25.sh 9.3: animal that used to decorate the Unix room
cat $IN93 | cut -c 1-2 | tr -d '\n'

#26.sh 9.4: four corners with E centered, for an "X" configuration
cat $IN94 | tr ' ' '\n' | grep "\"" | sed 4d | cut -d "\"" -f 2 | tr -d '\n'

#27.sh # 9.5: backwards running clock, in a backwards poem
cat $IN95 | rev | cut -c -3 | rev | grep -o '[A-Z][A-Z]*' | rev | tr -d '\n' | rev

#28.sh 9.6: Follow the directions for grep
cat $IN96 | tr ' ' '\n' | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'

#29.sh 9.7: Four corners
cat $IN97 | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'

#30.sh 9.8: TELE-communications
cat $IN98 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 2d | sed 3d | sed 4d | tr -c '[A-Z]' '\n' | tr -d '\n'

#31.sh 9.9:
cat $IN99 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 1d | sed 2d | sed 3d | sed 5d | tr -c '[A-Z]' '\n' | tr -d '\n'

#32.sh 10.1: count Turing award recipients while working at Bell Labs
cat $IN10 | sed 1d | grep 'Bell' | cut -f 2 | wc -l

#33.sh 10.2: list Turing award recipients while working at Bell Labs
cat $IN10 | sed 1d | grep 'Bell' | cut -f 2

#34.sh 10.3: extract Ritchie's username
cat $IN10 | grep 'Bell' | cut -f 2 | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

#35.sh 11.1: year Ritchie and Thompson receive the Hamming medal
cat $IN11 | grep 'UNIX' | cut -f 1

#36.sh 11.2: most repeated first name in the list?
cat $IN11 | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d
