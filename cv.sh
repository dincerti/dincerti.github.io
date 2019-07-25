#!/bin/sh

# pdf
/Applications/RStudio.app/Contents/MacOS/pandoc/pandoc cv-pandoc.md --template latex-cv-template.tex \
--variable author="Devin Incerti" --variable title="CV" --variable linkcolor="blue" \
--variable fontfamily="mathpazo" --variable fontfamilyoptions="sc, osf" --variable fontsize="11pt" \
--variable geometry="margin=1in" --variable jobtitle="Senior Data Scientist" \
--variable email="devin.incerti@gmail.com" --variable fontawesome="yes" \
--variable address="350 DNA Way, South San Francisco, CA 94080" --variable urlcolor="blue" \
--variable web="devinincerti.com" \
--pdf-engine /Library/TeX/texbin/pdflatex -o dincerti-cv.pdf

# html
/Applications/RStudio.app/Contents/MacOS/pandoc/pandoc cv-pandoc.md -o cv.html
echo '---\nlayout: default\n--- \n<p> In addition to the HTML version of my CV, 
a PDF version is available <a href="dincerti-cv.pdf"> here</a>.' > temp_file.csv
cat cv.html >> temp_file.csv
mv temp_file.csv cv.html

