#!/bin/bash
#SBATCH --mem=24000 # max is 64G on common partition, 64-240G on common-large
#SBATCH --output=/hpc/home/%u/FS6ADNI.%j.out
#SBATCH --error=/hpc/home/%u/FS6ADNI.%j.out

sub=$1  # cross-sectional step: RID_DATE, e.g. 60_20070709; longitudinal step: RID only, e.g. 60
step=$2 # "cross" for standard/cross-sectional step, "long" for template step (takes a while) + longitudinal step (shorter). must be run in that order!
		# "long" step assumes that cross step is finished for each of the subject's timepoints, and stored in $bidsDir/${sub}_* (RID_DATE name format, where DATE=YYYYMMDD)
img=$3  # full path to either single nifti or first dcm in the downloaded data folder ("cross" step only)

bidsDir=$H/Studies/ADNI/BIDS/derivatives/freesurfer_v6.0 # final output directory

# these only get used for downloaded folders that only contain one dicom and need to be converted before being passed to freesurfer
rawDataPath=$H/Studies/ADNI/raw
convertedDataPath=$H/Studies/ADNI/raw/dcm2nii

export SUBJECTS_DIR=/work/long/freesurfer_v6.0_ADNI # temporary working dir
export FREESURFER_HOME=/cifs/hariri-long/Scripts/Tools/FreeSurfer/freesurfer
export PATH=$FREESURFER_HOME/bin:$FREESURFER_HOME:$PATH


# QC requires a conda environment called freesurferQA with run_fsqc installed
condaPath=/hpc/home/ark19/miniconda3/etc/profile.d

echo "----JOB [$SLURM_JOB_ID] START [`date`] on HOST [$HOSTNAME]----"
echo "----CALL: $0 $@----"

module load compatbin # added 1/7/24 to fix OS compatibility issues after cluster update, with help from Tom Milledge 

for fsavg in fsaverage fsaverage_sym fsaverage3 fsaverage4 fsaverage5 fsaverage6; do
	if [[ ! -d $SUBJECTS_DIR/$fsavg ]]; then cp -r ${FREESURFER_HOME}/subjects/$fsavg $SUBJECTS_DIR; fi
done

if [[ $step == "cross" ]]; then

	if [[ ! -e $SUBJECTS_DIR/$sub/scripts/recon-all.done ]] || [[ -e $SUBJECTS_DIR/$sub/scripts/recon-all.error ]]; then

		rm -r $SUBJECTS_DIR/$sub/
		
		# some downloaded dirs only have a single dicom that seems to be a format that FS can't handle, but works to pre-convert to nii then pass to FS
		if [[ $(ls $(dirname $img)/*dcm | wc -l) -eq 1 ]]; then
			outdir=$(dirname $img | sed "s|$rawDataPath|$convertedDataPath|")
			mkdir -p $outdir
			dcm2niix -o $outdir $(dirname $img)
			img=$(ls $outdir/*gz | head -1)
		fi

		echo "RUNNING: $FREESURFER_HOME/bin/recon-all_noLink -subjid $sub -all -openmp 1 $str"
		$FREESURFER_HOME/bin/recon-all_noLink -subjid $sub -all -openmp 1 -i $img

	fi
	
	# run QC
	if [[ -e $SUBJECTS_DIR/$sub/scripts/recon-all.done ]] && [[ ! -e $SUBJECTS_DIR/$sub/scripts/recon-all.error ]]; then
		source $condaPath/conda.sh
		conda activate freesurferQA
		run_fsqc --subjects_dir $SUBJECTS_DIR --output_dir $SUBJECTS_DIR/$sub/QA --subjects $sub --screenshots --surfaces --skullstrip
	fi
	
	outdir=$SUBJECTS_DIR/$sub # for slurm output file
	
fi # end if step==cross

if [[ $step == "long" ]]; then

	# Check if template has been run, run if not
	if [[ -e $bidsDir/sub-${sub}_template/scripts/recon-all.done ]] && [[ ! -e $bidsDir/sub-${sub}_template/scripts/recon-all.error ]]; then 
		echo "!!!!!!!!!!!!!!!!!!!!!!!!! Completed template step already stored in $bidsdir! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		exit
	else
		if [[ -e $SUBJECTS_DIR/sub-${sub}_template/scripts/recon-all.done ]] && [[ ! -e $SUBJECTS_DIR/sub-${sub}_template/scripts/recon-all.error ]]; then 
			echo "!!!!!!!!!!!!!!!!!!!!!!!!! Completed template step already stored in $SUBJECTS_DIR! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			# mv $SUBJECTS_DIR/sub-${sub}_template $bidsDir
		else
		
			rm -rf $SUBJECTS_DIR/sub-${sub}_template

			##Set up directories
			str=""
			for f in $(ls -d $bidsDir/${sub}_*); do 
				if [[ ! -e $f/scripts/recon-all.done ]] || [[ -e $f/scripts/recon-all.error ]]; then
					echo "!!! Cross sectional step for $f not completed, exiting !!!"
					exit
				else
					scanID=$(basename $f)
					rm -rf $SUBJECTS_DIR/$scanID
					mkdir -p $SUBJECTS_DIR/$scanID
					cp -r $f/* $SUBJECTS_DIR/$scanID/
					str="$str -tp $scanID"
					echo "ADDING TIMEPOINT: $scanID"
				fi
			done
			
			# Create an unbiased template from all time points for each subject and process it with recon-all:
			echo "RUNNING: ${FREESURFER_HOME}/bin/recon-all_noLink -base sub-${sub}_template $str -all"
			${FREESURFER_HOME}/bin/recon-all_noLink -base sub-${sub}_template $str -all
			
		fi # end if not already completed on /work
	fi # end if not already completed on H:
			
	# Now run Step 3: "-long" longitudinally process all timepoints:
	for f in $(ls -d $bidsDir/${sub}_*); do
	
		scanID=$(basename $f)
		
		dir=$bidsDir/$scanID.long.sub-${sub}_template
		if [[ -e $dir/scripts/recon-all.done ]] && [[ ! -e $dir/scripts/recon-all.error ]]; then
			echo "!!!!!!!!!!!!!!!!!!!!!!!!! Completed longitudinal step for $scanID already stored in $bidsDir! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		else
			dir=$SUBJECTS_DIR/$scanID.long.sub-${sub}_template
			if [[ -e $dir/scripts/recon-all.done ]] && [[ ! -e $dir/scripts/recon-all.error ]]; then
				echo "!!!!!!!!!!!!!!!!!!!!!!!!! Completed longitudinal step for $scanID already stored in $SUBJECTS_DIR! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			else		
				
				rm -r $dir
				${FREESURFER_HOME}/bin/recon-all_noLink -long $scanID sub-${sub}_template -all	
				
			fi # end if not already completed in temp dir
		fi # end if not already completed in final dir
		
	done
	
	outdir=$SUBJECTS_DIR/sub-${sub}_template # for slurm output file
	
fi # end if step==long

# -- BEGIN POST-USER --
echo "----JOB [$SLURM_JOB_ID] STOP [`date`]----"
mv /hpc/home/$USER/FS6ADNI.$SLURM_JOB_ID.out $outdir/FS6ADNI.$SLURM_JOB_ID.out
# -- END POST-USER --
