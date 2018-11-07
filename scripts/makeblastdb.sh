#!/usr/bin/env bash
echo $0
echo $1
makeblastdb -in $1 -dbtype prot