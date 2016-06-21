rem ffmpeg -y -f image2 -r 15 -i sol_%%03d.png -c:v libx264 -r 25 sol.mp4
rem ffmpeg -y -f image2 -r 15 -i sol_%%05d.png  -c:v mpeg4 -qscale:v 1 -r 15 sol.avi
rem ffmpeg -y -f image2 -r 15 -i solz_%%03d.png  -r 15 solz.avi
rem ffmpeg -y -f image2 -r 15 -i phase1_%%03d.png  -r 15 phase.wmv
rem ffmpeg -y -f image2 -r 15 -i phase_%%05d.png  -c:v mpeg4 -qscale:v 1 -r 15 phase.avi
rem ffmpeg -y -f image2 -r 15 -i phasez_%%03d.png -r 15 phasez.avi
rem ffmpeg -y -f image2 -r 15 -i phase_%%03d.png -c:v libx264 -r 25 phase.mp4
rem ffmpeg -y -f image2 -r 15 -i phase%%d.png  -r 15 phase.avi
rem ffmpeg -y -f image2 -r 15 -i core_%%03d.png -c:v mjpeg -r 15 -qscale 2 core.avi
rem ffmpeg -y -f image2 -r 15 -i core_%%05d.png  -c:v mpeg4 -qscale:v 1 -r 15 core.avi
ffmpeg -y -f image2 -r 15 -i combined_%%05d.png  -c:v mpeg4 -qscale:v 1 -r 15 combined.avi