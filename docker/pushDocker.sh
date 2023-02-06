#!/bin/bash
#
docker build --progress=plain -t fcunial/lr-genotyping .
docker push fcunial/lr-genotyping
