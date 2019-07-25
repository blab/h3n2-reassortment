rm *.aux
rm *.log
rm *.out
rm *.pdf

pdflatex response1.tex
pdflatex response1.tex

open response1.pdf
