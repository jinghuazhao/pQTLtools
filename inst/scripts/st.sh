#!/usr/bin/bash

bash Olink_NGS.sh
ls | parallel -C' ' 'R --no-save <{}'
mv *rda ../../data
