rm -rf ~/temp/INLAspacetime*
mkdir ~/temp/INLAspacetime/
cp -r DESCRIPTION NAMESPACE R/ man/ src/ demo/ vignettes/ ~/temp/INLAspacetime/
 cd vignettes/
 mkdir ~/temp/INLAspacetime/vignettes/
 cp -r figure preamble.tex references.bib *.Rmd ~/temp/INLAspacetime/vignettes/
cd ~/temp/
rm INLAspacetime/R/cgeneric_ast.R
rm INLAspacetime/man/cgeneric_ast.Rd
rm INLAspacetime/R/cgeneric_nngp.R
rm INLAspacetime/man/cgeneric_nngp.Rd
R CMD build INLAspacetime
R CMD check INLAspacetime_*.tar.gz --as-cran


