
for in1 in *A\ horizontal.jpg; do
	in2=$(echo "$in1" | sed 's/A horizontal.jpg/B horizontal.jpg/g')
	out=$(echo "$in1" | sed 's/A horizontal.jpg/combined.jpg/g')
	magick "$in1" "$in2" +append "$out"
done

outdir=combined
mkdir $outdir 1>/dev/null
mv *\ combined.jpg ${outdir}