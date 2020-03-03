#!/usr/bin/env bash
# Sample command line for imghist

if [ "$1" == "-c" ]; then
  ./sqrs -r test-100
  ./lzw -c test-100.qrs
#  ./lzw -d test-100.qrs.lzw
elif [ "$1" == "-e" ]; then
  ./sqrs -r test-100
  ./aes -e 128 -i test-100.qrs -o enc-test-100
#  ./aes -d 128 -i enc-test-100 -o test-100.qrs
elif [ "$1" == "-ce" ]; then
  ./sqrs -r test-100
  ./lzw -c test-100.qrs
  ./aes -e 128 -i test-100.qrs.lzw -o enc-test-100-lzw
#  ./aes -d 128 -i enc-test-100-lzw -o test-100.qrs.lzw
#  ./lzw -d test-100.qrs.lzw
else
  ./sqrs -r test-100
fi
