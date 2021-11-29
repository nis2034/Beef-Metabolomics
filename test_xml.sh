#!/bin/bash
acc=($(grep -oP '(?<=email>)[^<]+' "/my_path/spam.xml"))

for i in ${!emails[*]}
do
  echo "$i" "${emails[$i]}"
  # instead of echo use the values to send emails, etc
done

