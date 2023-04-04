#!/bin/bash
################################################
# Start a self-resubmitting simulation.
################################################

# ID number for run
CASE_NAME="Couple"
JOBNO=319

# clean run directory and link all required files
./Coupled_prepare_run.sh

# record start times
#TIMEQSTART="$(date +%s)"
#echo Start-time `date` >> ../run/times

#ENV variables

totalruntime=3888000000 # 125 years long run total time of run: build into script
coupleruntime=1728000  #20 days time period between couplings no melt
ncouplesteps=$((totalruntime/coupleruntime)) #number of coupling runs

WAVI_deltaT_secs=1728000 #WAVI timestep in seconds (note that this will be converted to year in WAVI 3600*24*365)
WAVI_couplerunsteps=$((coupleruntime/WAVI_deltaT_secs)) #number of MITgcm time steps in a couple run

MITgcm_deltaT=200 #MITgcm timestep
MITgcm_couplerunsteps=$((coupleruntime/MITgcm_deltaT)) #number of MITgcm time steps in a couple run

#output frequency for models
WAVI_output=1728000  #note needs to be multiple of timestep
MITgcm_output=1728000  #only 7 outputs

couple_step=1 #set to 1 if new run
Ice_ran=0 #If Oce model has ran in loop or not

# Directory to save copies of pickups in
OUT_DIR="pickups"

Run_dir="/data/hpcdata/users/davbet33/mitgcm/cases/$CASE_NAME$JOBNO/run"


# submit the job chain
sbatch --export=CASE_NAME=$CASE_NAME,JOBNO=$JOBNO,totalruntime=$totalruntime,coupleruntime=$coupleruntime,ncouplesteps=$ncouplesteps,WAVI_deltaT_secs=$WAVI_deltaT_secs,WAVI_couplerunsteps=$WAVI_couplerunsteps,MITgcm_deltaT=$MITgcm_deltaT,MITgcm_couplerunsteps=$MITgcm_couplerunsteps,,WAVI_output=$WAVI_output,MITgcm_output=$MITgcm_output,couple_step=$couple_step,Ice_ran=$Ice_ran,OUT_DIR=$OUT_DIR,Run_dir=$Run_dir,HECACC=$HECACC,TIMEQSTART=$TIMEQSTART,IMGNAME=$IMGNAME,IMGNAME2=$IMGNAME2,JDEPOT=$JDEPOT,JOBNO=$JOBNO,TIMEQSTART=$TIMEQSTART  -J Coup_$JOBNO -A$HECACC ./Coupled_run.sh
