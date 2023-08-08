#!/bin/bash
# plots the output of Daniel Trugman's stress inversion outputs
# within a selected time range
# using GMT6
# last updated 3 Aug 2023
#
#                   
combined_kfile=${1-stress_inversion_binX-0.50_strideX-0.50_binY-0.50_strideY-0.50_binT-0.50_strideT-0.50.f32.txt}

tmin=${2-2010}                  # time range
tmax=${3-2011}

reg='-R129/145/31/44'
proj=-JM7	
msize=0.3

ostring=${combined_kfile%.*}	
ofile=$ostring.ps		# output file name, PS, temporary
ofile_pdf=$ostring.pdf	

gawk -v scale="$scale" -v tmin="$tmin" -v tmax="$tmax" 'FNR>1 && ($3>=tmin) && ($3<tmax)' "$combined_kfile" > "$tmpn.kfile_dh"

# lon lat time count friction misfit shape s1a s1p s2a s2p s3a s3p s11 s12 s13 s22 s23 s33
gawk '{print($14,$17,$19,$15,$16,$18,$15,$16,$18)}' $tmpn.kfile_dh | gawk -f normrow.awk > $tmpn.norm
rm tmpn.norm.v
cp $tmpn.norm ./tmpn.norm.v

paste $tmpn.kfile_dh $tmpn.norm | \
            gawk '{s=sqrt(0.2/(1e-3 + 0.2));print($1,$2,($17+$19)/2/$20,$14*s,$17*s,$19*s,$15*s,$16*s,$18*s,22,$1,$2)}' > $tmpn.dat

gawk '{print($8)}' $tmpn.kfile_dh > $tmpn.azi
gawk '{print($7)}' $tmpn.kfile_dh > $tmpn.R
gawk '{print($1,$2)}' $tmpn.kfile_dh > $tmpn.lonlat

paste $tmpn.lonlat $tmpn.R $tmpn.azi > $tmpn.dat2

gmt pscoast $reg $proj -Dl+f -W.25 -K -P -S100 -G200 > $ofile
gmt makecpt -T0/1/0.1 -Croma > $tmpn.cpt
gmt makecpt -T-0.4/0.4/0.05 -D -Croma > $tmpn.def.cpt
# gmt psmeca $tmpn.dat $reg $proj  -L0.25 -W0.25 -D-1e5/1e5 -Z$tmpn.def.cpt \
# 	       -Sm$msize+l+s1e22 -O -K  >> $ofile

gmt psmeca tmpn.dat $reg $proj  -L0.25 -W0.25 -D-1e5/1e5 -C$tmpn.def.cpt \
               -Sm$msize+l+s1e22 -O -K  >> $ofile

gawk '{print($1,$2,$3,$4,0.3)}' $tmpn.dat2 | \
   gmt psxy -C$tmpn.cpt -SVB0.04/0.01/0.01 $reg $proj -Gblack -O -K >> $ofile

gmt psscale -D1.5/4.5/3/.2 -C$tmpn.cpt -B.2/:"S": -O -K >> $ofile
#gmt psscale -D2.1/4.5/3/.2 -O -K -C$tmpn.def.cpt -B.2/:"@~s@~@-m@-": >> $ofile

gmt psbasemap $reg $proj -O -Ba5f1WeSn >> $ofile # add labels and end plot
gmt psconvert -Tf -A+m0.1 $ofile		      # convert to PDF and fix BB
echo $0: output in $ostring.pdf
rm $ofile
rm $tmpn
