#!/usr/bin/env bash
# Sample command line for imghist

./aes -e 128 -i README -o encREADME
./aes -d 128 -i encREADME -o decREADME