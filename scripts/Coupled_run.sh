#!/bin/bash
################################################################################
# Run the model for as long as we can, then prepare for a restart and submit the next job.
################################################################################
#
#SBATCH --chdir=../run/
#SBATCH --output=job.out
##SBATCH --mail-type=begin,end,fail,requeue    ## When to email you
##SBATCH --mail-user=davbet33@bas.ac.uk        ## Your email address
#SBATCH --time=24:00:00                       ## Maximum CPU time
#SBATCH --nodes=1                             ## Run on 1 node
#SBATCH --ntasks-per-node=126                  ## Run 44 tasks per node
#SBATCH --cpus-per-task=1                ## Run 44 tasks per node
##SBATCH --reservation=shortqos
#SBATCH --partition=standard                    ## Which Partition/Queue to use
#SBATCH --qos=standard                    ## Which Partition/Queue to use
##SBATCH --qos=short                    ## Which Partition/Queue to use

# function to return number of seconds left in this job
# the squeue command returns either hh:mm:ss or mm:ss
# so handle both cases.
# We should add in 1-00:00:00 for a day

# pass the julia depot path through singularity
export SINGULARITYENV_JULIA_DEPOT_PATH="/opt/julia"


#source /etc/profile.d/modules.sh
#module purge
echo "Loading Modules"
#Modules for MITgcm

###module load epcc-job-env

#module load hpc/hdf5/gcc/1.10.4
#module load hpc/mvapich2/gcc/2.2
#module load hpc/netcdf/gcc/4.4.1.1
#module load hpc/gcc/7.2.0
#module load hpc/openmpi/gcc

#Modules for WAVI

#module load hpc/python/conda-python-3.7.3


# function to return number of seconds left in this job
# the squeue command returns either hh:mm:ss or mm:ss
# so handle both cases.
# We should add in 1-00:00:00 for a day


function hmsleft()
{
        local lhms
        lhms=$(squeue  -j $SLURM_JOB_ID -O TimeLeft | tail -1)
        echo $lhms
}
function secsleft() {
    if [[ ${#hms} < 6 ]]
    then
        echo secs=$(echo $hms|awk -F: '{print ($1 * 60) + $2 }')
    else
        echo secs=$(echo $hms|awk -F: '{print ($1 * 3600) + ($2 * 60) + $3 }')
    fi
}



# start global timer
timeqend="$(date +%s)"
timestart="$(date +%s)"
echo >> times
echo Run start `date` >> times
hms=$(hmsleft)
echo Walltime left is $hms>>walltime
rem_secs=$(secsleft)  # function above
echo Walltime left in seconds is $rem_secs >> walltime
# Subtract 15 minutes
RUNTIME="$(($rem_secs-900))"
echo Will run for $RUNTIME sec >> walltime



echo "couple_step=", $couple_step
echo "dir", $Run_dir

echo "received from SLURM"  HECACC=$HECACC,IMGNAME=$IMGNAME2,JDEPOT=$JDEPOT,JOBNO=$JOBNO,TIMEQSTART=$TIMEQSTART

edit11="Couple_run_end = $totalruntime"
sed "s/Couple_run_end = unset/$edit11/" driver_updated.jl > driver.temp6 
mv driver.temp6 driver_updated.jl


#if [ $couple_step == 1 ] && [ $Oce_ran == 0 ]; then      
if [ $couple_step == 1 ]; then      

# Create output directory if it doesn't exist (to save copies of the pickup files)
  if [ ! -d $OUT_DIR ] ; then
    mkdir $OUT_DIR
  fi




echo "edit files"
edit=" ntimesteps=$MITgcm_couplerunsteps"
edit2=" streamice_couple_time= $coupleruntime",
edit3=" deltaT= $MITgcm_deltaT",
edit4=" dt_s = $WAVI_deltaT_secs"
edit5="end_time_s=$coupleruntime"
edit6=" pChkptFreq = $coupleruntime"
edit7=" WAVI_run_dir = $Run_dir"
edit8="chkpt_freq_s = $coupleruntime"
edit9="pchkpt_freq_s = $coupleruntime"
edit10="output_freq_s = $WAVI_output"

edit12=" frequency(1) = $MITgcm_output.,"
edit13=" frequency(2) = $MITgcm_output.,"
edit14=" frequency(3) = $MITgcm_output.,"
edit15=" frequency(4) = $MITgcm_output.,"
edit16=" frequency(5) = $MITgcm_output.,"
edit17=" frequency(6) = $MITgcm_output.,"
edit18=" frequency(7) = $MITgcm_output.,"
edit19=" frequency(8) = $MITgcm_output.,"

edit20=" dt_coup = $MITgcm_deltaT"

sed "s/.*ntimesteps.*/$edit/" data > data.temp1 
sed "s/.*streamice_couple_time.*/$edit2/" data.streamice > data.temp2
sed "s/.*deltaT.*/$edit3/" data.temp1 > data.temp3 
sed "s/.*pChkptFreq.*/$edit6/" data.temp3 > data.temp4 
sed "s/dt_s = unset/$edit4/" driver_updated.jl > driver.temp 
sed "s/end_time_s = unset/$edit5/" driver.temp > driver.temp2 
sed "s/chkpt_freq_s = unset/$edit8/" driver.temp2 > driver.temp3 
sed "s/pchkpt_freq_s = unset/$edit9/" driver.temp3 > driver.temp4 
sed "s/output_freq_s = unset/$edit10/" driver.temp4 > driver.temp5 
sed "s/dt_coup = unset/$edit20/" driver.temp5 > driver.temp7
sed "s/ frequency(1) = unset.,/$edit12/" data.diagnostics > data.temp10
sed "s/ frequency(2) = unset.,/$edit13/" data.temp10 > data.temp11
sed "s/ frequency(3) = unset.,/$edit14/" data.temp11 > data.temp12
sed "s/ frequency(4) = unset.,/$edit15/" data.temp12 > data.temp13
sed "s/ frequency(5) = unset.,/$edit16/" data.temp13 > data.temp14
sed "s/ frequency(6) = unset.,/$edit17/" data.temp14 > data.temp15
sed "s/ frequency(7) = unset.,/$edit18/" data.temp15 > data.temp16
sed "s/ frequency(8) = unset.,/$edit19/" data.temp16 > data.temp17

mv data.temp4 data
mv data.temp2 data.streamice
mv data.temp17 data.diagnostics
mv driver.temp7 driver_updated.jl

#Model timers 

WAVI_elapsedtotal=0
MITgcm_elapsedtotal=0



fi

while [ $couple_step -le $ncouplesteps ]; do

 echo "COUPLED STEP $couple_step"


if [ $Ice_ran == 0 ]; then

 
#SETUP and run WAVI

# function to return number of seconds left in this job
function hmsleft()
{
        local lhms
        lhms=$(squeue  -j $SLURM_JOB_ID -O TimeLeft | tail -1)
        echo $lhms
}
function secsleft() {
    if [[ ${#hms} < 6 ]]
    then
        echo secs=$(echo $hms|awk -F: '{print ($1 * 60) + $2 }')
    else
        echo secs=$(echo $hms|awk -F: '{print ($1 * 3600) + ($2 * 60) + $3 }')
    fi
}

hms=$(hmsleft)
echo Walltime left is $hms>>walltime
rem_secs=$(secsleft)  # function above
echo Walltime left in seconds is $rem_secs >> walltime
# Subtract 15 minutes
RUNTIME="$(($rem_secs-900))"
echo Will run for $RUNTIME sec >> walltime





 # start WAVI timer
 WAVI_timestart="$(date +%s)"


 echo "Running WAVI"

 timeout $RUNTIME singularity exec -B ${JDEPOT}:/opt/julia,/mnt/lustre/a2fs-work2/work/n02/n02/davbet33/wavi/WAVIhpc/cases/,/work/n02/n02/davbet33/wavi/WAVIhpc/cases/ ${IMGNAME2} julia driver_updated.jl

# Get the exit code
  OUT=$?
  echo $OUT

 echo "WAVI ran"



 # end WAVI timer
 WAVI_timeend="$(date +%s)"
 WAVI_elapsedrun="$(expr $WAVI_timeend - $WAVI_timestart)"
 WAVI_elapsedtotal="$(expr $WAVI_elapsedtotal + $WAVI_elapsedrun)"


 for file in PChkpt_*.jld2; do
    [[ $file -nt $PICKUP_FILE ]] && PICKUP_FILE=$file
 done
  # Extract the middle bit of this filename
  PICKUP=${PICKUP_FILE#PChkpt_}
  PICKUP=${PICKUP%.jld2}

  re='^[0-9]+$'
  if [[ $PICKUP =~ $re ]]; then
# Save the timestep, with any leading zeros removed
  NITER0=$(echo $PICKUP | sed 's/^0*//')


  #edit the driver namelist: replace the first instance of niter0 = * with appropriate checkpoint number
  NITER0_LINE="niter0 = $NITER0"
  echo $NITER0_LINE
  sed -i '0,/.*niter0.*/s//'"$NITER0_LINE"'/' driver_updated.jl


fi




#if time ran out on Ice model then break loop
if [ $OUT == 124 ] || [ $OUT == 125 ]; then

  echo 'WAVI ran out of time'

  break

fi


  couple_plus=$((couple_step+1))
  Ice_end=$((coupleruntime * couple_plus))
  edit14="end_time_s=$Ice_end"
  sed -i '0,/.*end_time_s.*/s//'"$edit14"'/' driver_updated.jl

 #Set Ice switch
     Ice_ran=1



fi
#######################################################################

 #Run MITgcm
 # start MITgcm timer
 MITgcm_timestart="$(date +%s)"

 echo "Running MITgcm"
# timeout $RUNTIME  mpirun -np $SLURM_NTASKS mitgcmuv
 timeout $RUNTIME srun --distribution=block:block --hint=nomultithread ./mitgcmuv

# Get the exit code
 OUT=$?
 echo $OUT
 echo "MITgcm ran"

# end MITgcm timer
 MITgcm_timeend="$(date +%s)"
 MITgcm_elapsedrun="$(expr $MITgcm_timeend - $MITgcm_timestart)"
 MITgcm_elapsedtotal="$(expr $MITgcm_elapsedtotal + $MITgcm_elapsedrun)"



#if time ran out on Ocean model then break loop
if [ $OUT == 124 ] || [ $OUT == 125 ]; then

  echo 'MITgcm ran out of time'
  break

fi

#define nitero variable from looking at pickup files?


  for file in pickup.*.data; do
    [[ $file -nt $PICKUP_FILE ]] && PICKUP_FILE=$file
  done
  # Extract the middle bit of this filename
  PICKUP=${PICKUP_FILE#pickup.}
  PICKUP=${PICKUP%.data}

  re='^[0-9]+$'

  NITER0=$(echo $PICKUP | sed 's/^0*//')

# Edit the "data" namelist
  # Update the line which sets niter0 and is uncommented
  sed -i '/^ niter0/c\ niter0='"$NITER0"',' data

 # Save copies of pickup files (ocean and sea ice, data and meta)
   cp pickup.$PICKUP.data pickup.$PICKUP.meta pickup_shelfice.$PICKUP.data pickup_shelfice.$PICKUP.meta pickup_streamice.$PICKUP.data pickup_streamice.$PICKUP.meta $OUT_DIR/


 
 Ice_ran=0
 couple_step=$((couple_step+1))

done

# end global timer
timeend="$(date +%s)"
elapsedtotal="$(expr $timeend - $timestart)"
echo >> times
echo Run end `date` >> times
echo Run-time seconds $elapsedtotal >> times
echo WAVI Run-time seconds $WAVI_elapsedtotal >> times
echo MITgcm Run-time seconds $MITgcm_elapsedtotal >> times

cp job.out job.out$((couple_step-1))


if [ $OUT == 0 ]; then

  # Simulation completed

  echo 'job chain: finished'

elif [ $OUT == 124 ] || [ $OUT == 125 ]; then

# submit the job chain
  echo 'submit next job'

sbatch --export=CASE_NAME=$CASE_NAME,JOBNO=$JOBNO,totalruntime=$totalruntime,coupleruntime=$coupleruntime,ncouplesteps=$ncouplesteps,WAVI_deltaT_secs=$WAVI_deltaT_secs,WAVI_couplerunsteps=$WAVI_couplerunsteps,MITgcm_deltaT=$MITgcm_deltaT,MITgcm_couplerunsteps=$MITgcm_couplerunsteps,,WAVI_output=$WAVI_output,MITgcm_output=$MITgcm_output,couple_step=$couple_step,Ice_ran=$Ice_ran,OUT_DIR=$OUT_DIR,Run_dir=$Run_dir,HECACC=$HECACC,TIMEQSTART=$TIMEQSTART,IMGNAME=$IMGNAME,IMGNAME2=$IMGNAME2,JDEPOT=$JDEPOT,JOBNO=$JOBNO,TIMEQSTART=$TIMEQSTART,elapsedtotal=$elapsedtotal,WAVI_elapsedtotal=$WAVI_elapsedtotal,MITgcm_elapsedtotal=$MITgcm_elapsedtotal    -J Coup_$JOBNO -A$HECACC ../scripts/Coupled_run.sh


fi
