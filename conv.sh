#avconv -y -f image2 -r 15 -start_number 1 -i sol_%05d.png   -c:v mpeg4 -qscale:v 10 -r 15 sol.avi
#avconv -y -f image2 -r 15 -start_number 1 -i phase_%05d.png   -c:v mpeg4 -qscale:v 10 -r 15 phase.avi
avconv -y -f image2 -r 15 -start_number 1 -i core_%05d.png   -c:v mpeg4 -qscale:v 5 -r 15 core.avi
#avconv -y -f image2 -r 15 -start_number 1 -i combined_%05d.png   -c:v mpeg4 -qscale:v 5 -r 15 combined.avi
#avconv -y -f image2 -r 15  -i vel_%05d.png   -c:v mpeg4 -qscale:v 10 -r 15 vel.avi

pdflatex tp.tex
#convert -density 150 tp.pdf -resize 1200x900! tp.png
#pdftoppm -png -scale-to-x 800 -scale-to-y 500 tp.pdf > tp.png
pdftoppm -png -scale-to 800 -y 60 -H 500 tp.pdf > tp.png

avconv -y -loop 1 -i tp.png -t 00:00:02 -c:v mpeg4 -qscale:v 5 tp.avi

avconv -y -i "concat:tp.avi|core.avi" -c copy combined_tp.avi

#rm tp.aux
#rm tp.log
#rm tp.pdf
#rm tp.png
#rm tp.avi
