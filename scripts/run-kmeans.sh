#!/usr/bin/env bash
# Sample command line for imghist

if [ "$1" == "-c" ]; then
  ./kmeans -i test-kmeans.txt -c 10 -o test-kmeans.out.txt
  ./lzw -c test-kmeans.out.txt
#  ./lzw -d test-kmeans.out.txt.lzw
elif [ "$1" == "-e" ]; then
  ./kmeans -i test-kmeans.txt -c 10 -o test-kmeans.out.txt
  ./aes -e 128 -i test-kmeans.out.txt -o enc-test-kmeans-out
#  ./aes -d 128 -i enc-test-kmeans-out -o test-kmeans.out.txt
elif [ "$1" == "-ce" ]; then
  ./kmeans -i test-kmeans.txt -c 10 -o test-kmeans.out.txt
  ./lzw -c test-kmeans.out.txt
  ./aes -e 128 -i test-kmeans.out.txt.lzw -o enc-test-kmeans-out-lzw
#  ./aes -d 128 -i enc-test-kmeans-out-lzw -o test-kmeans.out.txt.lzw
#  ./lzw -d test-kmeans.out.txt.lzw
else
  ./kmeans -i test-kmeans.txt -c 10 -o test-kmeans.out.txt
fi
