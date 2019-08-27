rm *.log
rm *.aux
rm *.out
rm *.bbl
rm *.bcf
rm *.blg
rm .DS_Store
rm h3n2_reassortment.pdf
osascript -e 'quit app "Preview"'

pdflatex h3n2_reassortment.tex
bibtex h3n2_reassortment
pdflatex h3n2_reassortment.tex
pdflatex h3n2_reassortment.tex
open h3n2_reassortment.pdf
