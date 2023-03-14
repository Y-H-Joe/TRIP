#!/bin/bash
SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
tail -n +2 ${SHELL_FOLDER}/*.txt | cut -f 7 | tr ";" "\t" | tee ${SHELL_FOLDER}/url.0 | cut -f 1 > ${SHELL_FOLDER}/url.1
cat ${SHELL_FOLDER}/url.0 | cut -s -f 2 > ${SHELL_FOLDER}/url.2
cat ${SHELL_FOLDER}/url.1 > ${SHELL_FOLDER}/url
cat ${SHELL_FOLDER}/url.2 >> ${SHELL_FOLDER}/url
rm ${SHELL_FOLDER}/url.1
rm ${SHELL_FOLDER}/url.2
rm ${SHELL_FOLDER}/url.0
