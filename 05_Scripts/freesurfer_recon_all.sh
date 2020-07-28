########################################################
#
# Freesurfer T1 segmentation
# Kayti Keith - 6/4/19
#
########################################################

export input=/path-to-data/

for ID in Subj1 Subj2 Subj3 ... Subjn ; do
recon-all -all -i $input/${ID}/T1.nii -sd $input/${ID} -subjid ${ID} 
done