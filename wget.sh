#!/bin/bash
SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
cat ${SHELL_FOLDER}/url_filtered | while read line
do
    wget --retry-connrefused -t 10  --timeout=300 -q -c $line --directory-prefix=${SHELL_FOLDER}
done
