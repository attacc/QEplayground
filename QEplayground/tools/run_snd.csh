#!/bin/csh

setenv QE       /home/attacc/SOFTWARE/qe-6.1/bin/pw.x
setenv P2Y      /home/attacc/SOFTWARE/yambo-devel/bin/p2y
setenv YAMBO    /home/attacc/SOFTWARE/yambo-devel/bin/yambo

setenv YSETUP   ../../yambo.in_setup
setenv YOPTICS  ../../yambo.in_optics

set scf_name="hBN.scf.in"
set nscf_name="hBN.nscf.in"

foreach file (*.scf.in_M*)
   set folder = `echo $file | awk '{print substr($0,length($0)-6,7)}'`
   if ( -d $folder) then 
     echo "Skipping folder $folder ..."
   endif
   set scf  = "${scf_name}_${folder}"
   set nscf = "${nscf_name}_${folder}"
   echo " DFT inputs:  $scf and  $nscf "
   mkdir $folder
   cd $folder
   $QE -inp ../$scf
   $QE -inp ../$nscf
   cd bn.save
   $P2Y -N -F data-file.xml
   $YAMBO -M -N -F $YSETUP
   $YAMBO -M -N -F $YOPTICS
   cd ..
   cd ..
end
