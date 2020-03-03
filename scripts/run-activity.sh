#!/usr/bin/env bash
# Sample command line for imghist

if [ "$1" == "-c" ]; then
  ./activity -r test-100.ts
  ./lzw -c test-activity-out.txt
#  ./lzw -d test-activity-out.txt.lzw
elif [ "$1" == "-e" ]; then
  ./activity -r test-100.ts
  ./aes -e 128 -i test-activity-out.txt -o enc-test-activity-out
#  ./aes -d 128 -i enc-test-activity-out -o test-activity-out.txt
elif [ "$1" == "-ce" ]; then
  ./activity -r test-100.ts
  ./lzw -c test-activity-out.txt
  ./aes -e 128 -i test-activity-out.txt.lzw -o enc-test-activity-out-lzw
#  ./aes -d 128 -i enc-test-activity-out-lzw -o test-activity-out.txt.lzw
#  ./lzw -d test-activity-out.txt.lzw
else
  ./activity -r test-100.ts
fi
