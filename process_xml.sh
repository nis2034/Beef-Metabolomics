#!/bin/bash
XML_FILE=split.aa
ACC=$(awk -F"[<>]" '/metabolite/ {f=1} f && /update_date/ {getline;print $3;f=0}' $XML_FILE)
DIR_PARENT=$(awk -F"[<>]" '/direct_parent/ {print $3}' $XML_FILE)
echo "$ACC"
echo "$DIR_PARENT"

