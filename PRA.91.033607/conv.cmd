rem d:\programs\ImageMagick-6.8.9-6\convert.exe tp.pdf -density 600x600 -quality 90 -resize 1200x900 file.png
rem ffmpeg -y -f image2 -r 15 -i combined_%%05d.png  -c:v mpeg4 -qscale:v 10 -r 15 combined.avi

latex tp.tex
dvipng -T 800pt,570pt tp.dvi
d:\programs\ImageMagick-6.8.9-6\convert.exe tp1.png -resize 1200x900! tp.png

ffmpeg -y -loop 1 -r 15 -i tp.png -t 00:00:02 -vcodec mpeg4 -an tp.avi

ffmpeg -y -i "concat:tp.avi|combined.avi" -c copy combined_tp.avi

del tp.aux
del tp.log
del tp.dvi
del tp1.png
rem del tp.png
del tp.avi