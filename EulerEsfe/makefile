all: 
	gcc EulerEsfe.c -o EulerEsfe -lgsl -lgslcblas -lm

gnuplot:
	true

plots:
	./plots.sh

video:
	ffmpeg -threads 3 -r 30 -pattern_type glob -i '*.png' -c:v libx264 out.mp4
	#ffmpeg -r 1 -pattern_type glob -i '*.png' -c:v libx264 out.mp4


clean: 
	# rm -f p data*.dat data*.png *.mp4
	rm -f p 
	find . -name "*.png" -print0 | xargs -0 rm
	find . -name "data*.dat" -print0 | xargs -0 rm
	find . -name "*.mp4" -print0 | xargs -0 rm

cleanpng:
	rm -f data*.png
