## This tutorial shows how to run voxel-based-lesion-to-symptom--mapping (VLSM) with Advanced Normalization Tools (ANTs and ANTsR)
# Requirements: Both ANTs and ANTsR (see https://github.com/stnava/ for how to install these)
# Test is done on Linux with the latest version of available on Oct 13 2014
#
# What is VLSM? It's a method trying to put in relation the presence or absence of lesion in a voxel with some behavioral measures.
# In practice, a lesion mask is created on each patient's brain damage (stroke, sclerosis, etc.). These masks are brought to a common template space. Thus, for each voxel we have a value of 0 or 1 for each subject, depending whether the subject was damaged or not. This creates a hypothetic population of two groups, one that is lesioned in that voxel, the other that is not. The simplest way to compare the two groups is to compare them about a behavioral measure using a t-test. We can also use linear regression, the result will be identical, but we can include more than one behavioral measure at once. We expect that, if the voxel is critical for the behavioral measure, the two groups will be different in the behavioral measure. Note that we are still talking about a single voxel. This rationale is repeated in each voxel.
# What is this example doing?
# There are five steps involved. However, if you want to save time and see how to do VLSM on a set of masks, go directly to STEP 5.
# STEP 1: take 20 brains from a publicly available dataset.
# STEP 2: draw some fictitious lesions or download them from github.
# STEP 3: calculate normalization of T1 images into template space.
# STEP 3: normalize T1 images into template
# STEP 4: apply transformation to lesion masks
# STEP 5: do a VLSM.

# your ANTs bin folder (antsRegistration and antsCorticalThickness.sh should be there)
MYANTSPATH =/path/to/your/ANTs/bin/

# define your data folder
mydata=/data/jag/username
cd $mydata


############################### STEP 1: get the data
wget http://mindboggle.info/data/mindboggle101/NKI-TRT-20_volumes.tar.gz
tar -zxvf NKI-TRT-20_volumes.tar.gz  # after this you have a folder called NKI-TRT-20_volumes with 20 subjects in it


############################### STEP 2: draw virtual lesions
# Open each subject's T1 and draw a lesion around the left inferior frontal cortex. Then save the file as tempseg.nii in the subject's folder.
# To draw lesions we can use ITKsnap or MRIcron, but remember, the lesion is a 3D object in multiple slices; don't expect to draw a single slice and call that a lesion.
# Alternatively you can get download some lesion masks I prepared for this example:
wget --no-check-certificate https://github.com/dorianps/VLSMwithANTs/blob/master/Lesions.tar.gz?raw=true
tar -zxvf Lesions.tar.gz

############################### STEP 3:  normalize T1 images into template
# This may be the most effortful step. On a single may run 6-40 hours per subject, but you can cut times by using antsRegistration instead of antsCorticalThickness.sh
# Also, the example runs all 20 subjects, but this is unnecesasry, only 10 will be used for normalization. (See below).
# All we need from this step is transformation files (Affine+Warp).
for i in {1..20}
do
   cd $mydata/NKI-TRT-20_volumes/NKI-TRT-20-$i/
	bash antsCorticalThickness.sh -d 3 -o OnOasis025 -w 0.25 \
	-a $mydata/NKI-TRT-20_volumes/NKI-TRT-20-$i/t1weighted.nii.gz \
	-e $mydata/template/T_template0.nii.gz \
	-f $mydata/template/T_template0_BrainCerebellumRegistrationMask.nii.gz \
	-p $mydata/template/Priors2/priors%02d.nii.gz \
	-m $mydata/template/T_template0_BrainCerebellumProbabilityMask.nii.gz \
	-t $mydata/template/T_template0_BrainCerebellum.nii.gz \
	| tee $mydata/NKI-TRT-20_volumes/NKI-TRT-20-$i/log.txt
done


############################### STEP 4: Apply transformations to lesion masks
# Before you start, create a folder called VLSM inside $mydata/NKI-TRT-20_volumes/
# We will switch to R for the remaining steps
# At this point we have calculated transformations from all subjects, will simply apply those.

R  # start R console
mydata=/data/jag/username  # define our base folder in R

# load ANTsR
library(ANTsR)

for (i in 1:10 ) {
	# only subjects 1-10 will be run, those for which we created virtual lesions
	# for subjects 11-20 we will create empty images later

	# switch working dir to this guy's folder
	setwd(paste(mydata,'/NKI-TRT-20_volumes/NKI-TRT-20-',i,'/', sep = ""))

	# read images
	lesion <- antsImageRead("tempseg.nii", 3)
	brainmask <- antsImageRead("OnOasis025BrainExtractionMask.nii.gz", 3)

	# your lesion was over the skull, lets take only the part inside the brain
	brainlesion <- antsImageClone(lesion)
	ImageMath(3, brainlesion, "m", lesion, brainmask)
	antsImageWrite(brainlesion, "brainlesion.nii.gz")


	# these are the transformation files, in the correct order, to bring subject to template
	mytx <- list()
	mytx$fwdtransforms[1] <- paste(mydata,"/NKI-TRT-20_volumes/NKI-TRT-20-",i, "/OnOasis025SubjectToTemplate1Warp.nii.gz", sep="")
	mytx$fwdtransforms[2] <- paste(mydata,"/NKI-TRT-20_volumes/NKI-TRT-20-",i, "/OnOasis025SubjectToTemplate0GenericAffine.mat", sep="")

	# apply transformations to the lesion mask
	brainlesiontemplate <- antsImageClone(brainlesion)  # this will be the lesion in template space
	template <- antsImageRead(paste(mydata,"/template/T_template0.nii.gz", sep=""), 3)  # this is simply the template T1
	brainlesiontemplate <- antsApplyTransforms(d="3", moving=brainlesion, fixed=template, t=mytx$fwdtransform, n="NearestNeighbor")
	ThresholdImage("3", brainlesiontemplate, brainlesiontemplate, 0.5, 1)  # binarize the the normalized lesion, only 0s and 1s

	# save the normalized lesion mask in the VLSM folder
	# we are saving files called subject_01.nii.gz, subject_02.nii.gz, etc. The leading zero is important to order the subjects automatically later
	antsImageWrite(brainlesiontemplate, paste(mydata,"/NKI-TRT-20_volumes/VLSM/subject_", sprintf("%02d", i), ".nii.gz", sep=""))

}

# now let's create empty images for subjects 11-20 (non lesioned group)
ThresholdImage("3", brainlesiontemplate, brainlesiontemplate, 2, 2)
for (i in 11:20 ) {
	antsImageWrite(brainlesiontemplate, paste(mydata,"/NKI-TRT-20_volumes/VLSM/subject_", sprintf("%02d", i), ".nii.gz", sep=""))

}

# Briefly exit R to create a VLSM mask: voxels where at least one subject was lesioned
q(save="no")
cd $mydata/NKI-TRT-20_volumes/VLSM/
AverageImages 3 LesionAverage.nii.gz 0 subject_*.nii.gz
ThresholdImage 3 LesionAverage.nii.gz vlsm_mask.nii.gz 0.001 1


############################### STEP 5: do a VLSM
# At this point you should have all subject_*.nii.gz files in the VLSM directory. If you don't have them, download them with the following commented lines:
# cd $mydata/NKI-TRT-20_volumes
# wget --no-check-certificate https://github.com/dorianps/VLSMwithANTs/blob/master/VLSM_folder.tar.gz?raw=true
# tar -zxvf VLSM_folder.tar.gz

R  # get back to R
mydata=/data/jag/username  # define our base folder in R
library(ANTsR)  # load ANTsR

setwd(paste(mydata,"/NKI-TRT-20_volumes/VLSM/", sep=""))  # switch to VLSM dir
imageList <- Sys.glob("s*nii.gz")  # get the filenames
mask <- antsImageRead("vlsm_mask.nii.gz",3)  # get the mask
mat <- imagesToMatrix(imageList,mask)  # form a matrix, rows=subjects, columns=voxels

# create 20 values 10+10
x <- rnorm(n=10, mean=10, sd=1)  # first 10 are lower, they're lesioned
y <- rnorm(n=10, mean=14, sd=1)  # second 10 are higher, no behavioral deficit
dx = c(x,y)
# create fake age with no difference between lesions
age <- rnorm( nrow(mat), mean=65, sd=15 )  # create some fake ages, without a pattern

# run VLSM, just a linear model with two variables
vlsm <- lm( mat ~ age + dx )  # a linear model for each voxel in "mat"
vlsm_interpretation <- bigLMStats( vlsm )  ## get info about the results

# look at the results
names( vlsm_interpretation )
age_var <- grep("age",rownames(vlsm_interpretation$beta.t))
dx_var <- grep("dx",rownames(vlsm_interpretation$beta.t))
agebetas <- vlsm_interpretation$beta.t[ age_var ,]
dxbetas <- vlsm_interpretation$beta.t[ dx_var,]
hist( agebetas )  # T values histogram distribution for age factor, no large T values in too many voxels
hist( dxbetas )  # T values histogram for our measure of interest, large negative T values in many voxels

# make an image and add betas to it
betas <- antsImageClone(mask)
betas[mask==1] <- -(vlsm_interpretation$beta.t[ age_var ,])  # we also invert negative to positive for easier view in ITKsnap or MRIcron
antsImageWrite(betas,"AgeBetas.nii.gz")
betas[mask==1] <- -(vlsm_interpretation$beta.t[ dx_var ,]) # we also invert negative to positive for easier view in ITKsnap or MRIcron
antsImageWrite(betas,"MyMeasureBetas.nii.gz")


# That's all. Open the template image and overlay "MyMeasureBetas.nii.gz" or "AgeBetas.nii.gz" to see the results.
