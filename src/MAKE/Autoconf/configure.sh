#!/bin/bash

autoconf
./configure
ls | grep -v configure.ac | grep -v *.sh | xargs rm -vR
