#!/bin/bash
for i in animation_frames/*.ppm;
  do name='echo "$i" | cut -d'.' -f1'
  echo "$name"
  convert "$i" "${name}.png"
done