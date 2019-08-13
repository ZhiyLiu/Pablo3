

#To use this file execute: sh execImage2Distance.sh <Filename binary image>
#example sh execImage2Distance.sh Right_Hippocampus/Right_Hippocampus.mhd

./Image2SignedDistanceMap ${1}

fn=${1%.*}

matlab_exec=/usr/local/MATLAB/R2010b/bin/matlab

X="path = pwd;
   path = strcat(path, '/MATLAB/');
   addpath(path);
   antiAliasWrapper('${fn}-ddm.mhd','${fn}',[0.5 0.5]);
   quit;"
echo ${X} > temp.m
cat temp.m
${matlab_exec} -nodisplay -r temp.m

rm temp.m 


