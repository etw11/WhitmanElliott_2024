# Step 2 for running ADNI ADSP data through freesurfer's longitudinal pipeline on a slurm cluster

SUBJECTS_DIR=/work/long/freesurfer_v6.0_ADNI             # temporary working dir
bidsDir=$H/Studies/ADNI/BIDS/derivatives/freesurfer_v6.0 # final output directory. longitudinal step will pull all scans for each subject from this dir
FSscript=$H/Studies/ADNI/Scripts/runFreesurfer6.0.sh     # full path to FS script

# move finished cross-sectional runs from temp dir to final bids dir
for f in `ls -d $SUBJECTS_DIR/[0-9]*`; do 
	if [ -e $f/scripts/recon-all.done ] && [ ! -e $f/scripts/recon-all.error ]; then 
		mv $f $bidsDir 2>/dev/null
	fi
done 

# run longitudinal step
cd $bidsDir 
for id in `ls -d [0-9]* | cut -d_ -f1 | sort | uniq`; do 
	sbatch $FSscript $id long
done 

