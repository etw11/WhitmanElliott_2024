#!/bin/bash

###Script to run longitudinal FS pipeline on UKB

ID=$1 ### Ex just 1809259 instead of 1809259_20263_2_0, the script will find the 2 FS sessions

outDir=/ncf/cnl07/Users/Max/UKB_longFS.v6.0

module load freesurfer/6.0.0-ncf
export SUBJECTS_DIR=$outDir
cd $SUBJECTS_DIR

tp1=$(ls -d ${outDir}/${ID}_20263_2_0)
tp2=$(ls -d ${outDir}/${ID}_20263_3_0)


####Make the template if it isn't already done
if [[ ! -f ${outDir}/${ID}_20263_template/stats/aseg.stats ]];then
    echo "Remaking template"
    rm -r ${outDir}/${ID}_20263_template
    recon-all -base ${ID}_20263_template -tp ${tp1} -tp ${tp2} -all
fi

###Then make the updated FS Long directories
jobSubmit1="module load freesurfer/6.0.0-ncf;export SUBJECTS_DIR=${outDir};cd $SUBJECTS_DIR;recon-all -long ${tp1} ${ID}_20263_template -all"
jobSubmit2="module load freesurfer/6.0.0-ncf;export SUBJECTS_DIR=${outDir};cd $SUBJECTS_DIR;recon-all -long ${tp2} ${ID}_20263_template -all"

###Learned the hard way that you have to have 001.mgz in the mri/orig folder otherwise the long pipeline won't run. I'm still confused about why it isn't already there...

if [[ ! -f ${outDir}/${ID}_20263_2_0.long.${ID}_20263_template/stats/aseg.stats ]];then
    echo "Remaking timepoint _2_0 longitudinal folder"
    rm -r ${outDir}/${ID}_20263_2_0.long.${ID}_20263_template
    mri_convert ${outDir}/${ID}_20263_2_0/mri/orig.mgz ${outDir}/${ID}_20263_2_0/mri/orig/001.mgz
    sbatch -p fasse -t 16:00:00 -o ~/swarmOutput/slurm.%j.out --mem=8G --wrap "$jobSubmit1"
fi

if [[ ! -f ${outDir}/${ID}_20263_3_0.long.${ID}_20263_template/stats/aseg.stats ]];then
    echo "Remaking timepoint _3_0 longitudinal folder"
    rm -r ${outDir}/${ID}_20263_3_0.long.${ID}_20263_template
    mri_convert ${outDir}/${ID}_20263_3_0/mri/orig.mgz ${outDir}/${ID}_20263_3_0/mri/orig/001.mgz
    sbatch -p fasse -t 16:00:00 -o ~/swarmOutput/slurm.%j.out --mem=8G --wrap "$jobSubmit2"
fi
