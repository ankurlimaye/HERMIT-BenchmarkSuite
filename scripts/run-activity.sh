#!/usr/bin/env bash
# Sample command line for imghist

if [ "$1" == "-c" ]; then
  ./activity -r test-100.ts
  ./lzw -c activity-out.txt
#  ./lzw -d activity-out.txt.lzw
elif [ "$1" == "-e" ]; then
  ./activity -r test-100.ts
  ./aes -e 128 -i activity-out.txt -o enc-activity-out
#  ./aes -d 128 -i enc-activity-out -o activity-out.txt
elif [ "$1" == "-ce" ]; then
  ./activity -r test-100.ts
  ./lzw -c activity-out.txt
  ./aes -e 128 -i activity-out.txt.lzw -o enc-activity-out-lzw
#  ./aes -d 128 -i enc-activity-out-lzw -o activity-out.txt.lzw
#  ./lzw -d activity-out.txt.lzw
else
  ./activity -r test-100.ts
fi