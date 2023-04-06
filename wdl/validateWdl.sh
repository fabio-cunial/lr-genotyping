#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l Regenotype.wdl -i inputs/private_Regenotype.json
java -jar ${WOMTOOL_PATH} validate -l Standardize.wdl -i inputs/private_Standardize.json
java -jar ${WOMTOOL_PATH} validate -l Sv2Igv.wdl
java -jar ${WOMTOOL_PATH} validate -l JediGraph.wdl
