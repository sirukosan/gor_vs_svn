#!/usr/bin/env bash
blastp -query $1 -db $2 -evalue 0.01 -out $3 -outfmt 6