#!/usr/bin/env bash
# Sample command line for imghist

if [ "$1" == "-c" ]; then
  ./imghist -r test-lena.bmp
  ./lzw -c out.bmp
#  ./lzw -d out.bmp.lzw
elif [ "$1" == "-e" ]; then
  ./imghist -r test-lena.bmp
  ./aes -e 128 -i out.bmp -o enc-out
#  ./aes -d 128 -i enc-out -o out.bmp
elif [ "$1" == "-ce" ]; then
  ./imghist -r test-lena.bmp
  ./lzw -c out.bmp
  ./aes -e 128 -i out.bmp.lzw -o enc-out
#  ./aes -d 128 -i enc-out -o out.bmp.lzw
#  ./lzw -c out.bmp.lzw
else
  ./imghist -r test-lena.bmp
fi