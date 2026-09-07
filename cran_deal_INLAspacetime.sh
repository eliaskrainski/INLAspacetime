rm -rf ~/temp/INLAspacetime*
mkdir ~/temp/INLAspacetime/
mkdir ~/temp/INLAspacetime/demo
mkdir ~/temp/INLAspacetime/man
mkdir ~/temp/INLAspacetime/R
mkdir ~/temp/INLAspacetime/src
mkdir ~/temp/INLAspacetime/vignettes/
cp DESCRIPTION NAMESPACE ~/temp/INLAspacetime/
cp demo/00Index demo/*R ~/temp/INLAspacetime/demo/
cp man/*Rd ~/temp/INLAspacetime/man/
cp R/*R ~/temp/INLAspacetime/R/
cp src/*h src/*c src/Makevars ~/temp/INLAspacetime/src/
cd vignettes
cp -r figures web preamble.tex references.bib *.Rmd ~/temp/INLAspacetime/vignettes/
cd ~/temp/
rm INLAspacetime/R/cgeneric_ast.R
rm INLAspacetime/man/cgeneric_ast.Rd
rm INLAspacetime/R/cgeneric_nngp.R
rm INLAspacetime/man/cgeneric_nngp.Rd
R CMD build INLAspacetime
R CMD check INLAspacetime_*.tar.gz --as-cran

