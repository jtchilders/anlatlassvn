#!/usr/bin/env bash

COMPARE_FILES_IN_THIS_FOLDER=$1
COMPARE_TO_THIS_FILE=$2

ls $COMPARE_FILES_IN_THIS_FOLDER/* | xargs -I file diff -i -s -q $COMPARE_TO_THIS_FILE file | grep identical

