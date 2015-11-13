avconv -y -f image2 -r 15 -start_number 0 -i sol_%05d.png   -c:v mpeg4 -qscale:v 1 -r 15 sol.avi
avconv -y -f image2 -r 15 -start_number 0 -i phase_%05d.png   -c:v mpeg4 -qscale:v 10 -r 15 phase.avi
avconv -y -f image2 -r 15 -start_number 0 -i combined_%05d.png   -c:v mpeg4 -qscale:v 10 -r 15 combined.avi
#avconv -y -f image2 -r 15 -i vel_%05d.png -c:v mpeg4 -qscale:v 10 -r 15 vel.avi
#avconv -y -f image2 -r 15 -i vel_%05d.png -s 600x450 -c:v h264 -preset veryslow -tune animation -crf 30 -r 15  vel.avi


#pdflatex tp.tex
#convert -density 150 tp.pdf -resize 1200x900! tp.png

#avconv -y -loop 1 -i tp.png -t 00:00:02 -vcodec mpeg4 tp.avi

#avconv -y -i "concat:tp.avi|combined.avi" -c copy combined_tp.avi

#rm tp.aux
#rm tp.log
#rm tp.pdf
#rm tp.png
#rm tp.avi