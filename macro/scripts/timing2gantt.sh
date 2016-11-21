#!/bin/bash

ostatus="ACTIVE"
ot="0.0"

sed -n '2,$p' < $1 | while IFS="," read macro_step macro_iter t status desc
do
  if [ "$ostatus" != "$status" ]
  then
    echo -e $2 '\t' $ot '\t' $t '\t' $ostatus
    ot=$t
    ostatus=$status
  fi
done
