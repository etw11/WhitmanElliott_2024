# Step 1 for running ADNI ADSP data through freesurfer's longitudinal pipeline on a slurm cluster

# full path to FS script
FSscript=$H/Studies/ADNI/Scripts/runFreesurfer6.0.sh

# this is the directory where all off the ADSP imaging data have been downloaded and unzipped.  
## Subdirectories should correspond to subjects, where the names are formatted as <site>_S_<RID>, e.g. 002_S_0295
ADSPdir=$H/Studies/ADNI/raw/ADSP-PHC.ALL    

# for each scan to run through FS, this file has 3 (tap or space separated) columns: 
## 1: <RID>_<DATE>, eg 295_20060418
## 2: ImageID, eg I45112 - this is the name of the bottom level dir containing dicoms/nii for this scan, and can be found in the ADSP csv file downloaded from IDA
## 3: full subject ID, formatted as <site>_S_<RID>, e.g. 002_S_0295 (can also be found in the csv download)	
IDfile=$H/Studies/ADNI/raw/ids_ADSPtoRun.txt 

# look up image paths and submit FS cross-sectional jobs
while read line; do 
	id=$(echo $line | awk '{print $1}') 
	img=$(echo $line | awk '{print $2}')
	s=$(echo $line | awk '{print $3}')
	img=$(ls $(find $ADSPdir/$s -name $img 2>/dev/null)/*{dcm,nii} | head -1)  
	sbatch $FSscript $id cross $img  
done < $IDfile 


