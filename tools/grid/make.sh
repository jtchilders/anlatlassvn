#!/usr/bin/env bash
# from: http://www.nikhef.nl/~janjust/proxy-verify/
gcc -o grid-proxy-verify grid-proxy-verify.c -lssl -lcrypto -std=c99

