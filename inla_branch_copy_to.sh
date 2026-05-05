cd src/
cp FUNCTIONS cgeneric_ar2ss.c cgeneric_ast.c cgeneric_barrier.c cgeneric_nngp.c cgeneric_sspde.c cgeneric_sstspde.c INLAspacetime.c INLAspacetime.h ../../inla_branch_INLAspacetime/src/
cd ../../inla_branch_INLAspacetime/
git add src/* 
git commit -m 'update inla branch'
git push
cd ../INLAspacetime/
