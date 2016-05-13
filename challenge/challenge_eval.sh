# Automation of validation and scoring
# Make sure you point to the directory where challenge.py belongs and a log directory must exist for the output
cd ./
#---------------------
#Validate submissions
#---------------------
python challenge.py -u SMC-RNA --send-messages --notifications validate --all >> log/score.log 2>&1

#--------------------
#Score submissions
#--------------------
python challenge.py -u SMC-RNA --send-messages --notifications score --all >> log/score.log 2>&1

#--------------------
#Cache submissions
#--------------------
#Fusion detection
python challenge.py -u SMC-RNA cache 5877348 syn6054684 
#Isoform Quantification
python challenge.py -u SMC-RNA cache 5952651 syn6054685 
