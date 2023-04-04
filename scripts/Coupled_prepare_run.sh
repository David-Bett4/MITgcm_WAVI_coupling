#!/bin/bash
################################################
# Clean out old results and link input files.
################################################

# Empty the run directory - but first make sure it exists!
if [ -d "../run" ]; then
  cd ../run
#  rm -rf *

  rm -rf *.bin
  rm -rf *.nc
  rm -rf *.jl
  rm -rf *.jl_old
  rm -rf *.zip
  rm -rf *.out
  rm -rf *.box
  rm -rf *.data
  rm -rf *.meta
  rm -rf data*
  rm -rf eedata
  rm -rf mitgcmuv
  rm -rf scratch*
  rm -rf *.jld
  rm -rf *.jld2
  rm -rf *.AS01*
  rm -rf *.m
  rm -rf *.mat
  rm -rf *.py
  rm -rf *.slurm
  rm -rf *.obw
  rm -rf *_mod4
  rm -rf *.init
  rm -rf *.jmd95z
  rm -rf pickups
  rm -rf job.out*
#  rm -rf julia-1.5.2
else
  echo 'There is no run directory'
  exit 1
fi

# Link everything from the input directory

ln -s ../input_WAVI/* . 
ln -s ../input_MITgcm/* . 

#copy files as editted in run directory

cp -f ../input_MITgcm/data .
cp -f ../input_MITgcm/data.streamice .
cp -f ../input_WAVI/driver_call_idealsiedsetup.jl .
cp -f ../input_WAVI/driver_combined.jl .
cp -f ../input_WAVI/driver_updated.jl .
cp -f ../input_WAVI/Run_test.jl .
cp -f ../input_WAVI/Run_test_cont.jl .
cp -f ../input_WAVI/Run_combined.jl .

# Link forcing files stored elsewhere

# Deep copy of the current version of wavi source directory

mkdir src
cp -rf $JDEPOT/packages/WAVI/. src/.


# Link executables
#ln -s ../../../julia-9d11f62bcb/bin/julia
#ln -s ../../../julia-9d11f62bcb/lib/lib* .
# ln -s ../../../julia-9d11f62bcb .
 ln -s ../build_MITgcm/mitgcmuv .
 ln -s ../scripts/mit2nc.py .
# ln -s /data/hpcdata/users/davbet33/mitgcm/MITgcm_c65z/utils/julia-1.5.2 .

#Link the image container
ln -s $IMGPATH .
ln -s $IMGPATH2 .


