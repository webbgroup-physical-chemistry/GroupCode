#! /bin/bash


for var in "$@" ; do
    i=$var
    echo $i

    ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 \
             -dSubsetFonts=true -dEmbedAllFonts=true \
             -dAutoFilterColorImages=false \
             -dAutoFilterGrayImages=false \
             -dColorImageFilter=/FlateEncode \
             -dGrayImageFilter=/FlateEncode \
             -dModoImageFilter=/FlateEncode \
			$i 
done
