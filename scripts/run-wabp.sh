#!/usr/bin/env bash
# Sample command line for imghist

if [ "$1" == "-c" ]; then
  ./wabp -r test-slp01a
  ./lzw -c test-slp01a.wabp
#  ./lzw -d test-slp01a.wabp.lzw
elif [ "$1" == "-e" ]; then
  ./wabp -r test-slp01a
  ./aes -e 128 -i test-slp01a.wabp -o enc-test-slp01a
#  ./aes -d 128 -i enc-test-slp01a -o test-slp01a.wabp
elif [ "$1" == "-ce" ]; then
  ./wabp -r test-slp01a
  ./lzw -c test-slp01a.wabp
  ./aes -e 128 -i test-slp01a.wabp.lzw -o enc-lzw-test-slp01a
#  ./aes -d 128 -i enc-lzw-test-slp01a -o test-slp01a.wabp.lzw
#  ./lzw -d test-slp01a.wabp.lzw
else
  ./wabp -r test-slp01a
fi