# $Id: limit_settings.py 661778 2015-04-20 22:23:21Z adye $
# Local settings for limit_submit

dsver       = opt.get("-V","1")
workspaces  = opt.get("-w","comb")

if grid:
    # ANALY_MPPMU,ANALY_LRZ are SLES clusters and are likely to have problems with pre-compiled stuff.
    # I found ANALY_GRIF-LPNHE,ANALY_MWT2,ANALY_SARA,ANALY_GRIF-LAL,ANALY_ARC,ANALY_GLASGOW,ANALY_DRESDEN,ANALY_ZA-UJ,ANALY_SFU_bugaboo had a rather high failure rate.
    # I found ANALY_HU_ATLAS_Tier2,ANALY_UCL,ANALY_INFN-GENOVA,ANALY_TRIUMF_CVMFS,ANALY_CAM,ANALY_GR-01-AUTH,ANALY_LPC particularly slow
    # ANALY_SLAC should be excluded if a Chirp server other than voatlas92.cern.ch is used. ANALY_SLAC has a hard-wired proxy connection to just voatlas92.
    # ANALY_TR-10-ULAKBIM failed lib job twice (17/6/12)
    # ANALY_NET2 is missing libkeyutils.so.1 (fixed 1/5/12, confirmed working 23/5/12)
    # ANALY_GRIF-IRFU has pilot error in build job (temp? failed on 1/5/12, worked the day before)
    excludeSite  = "ANALY_MPPMU,ANALY_LRZ,ANALY_GRIF-LPNHE,ANALY_MWT2,ANALY_SARA,ANALY_HU_ATLAS_Tier2"
    excludeSite += ",ANALY_UCL,ANALY_INFN-GENOVA,ANALY_GRIF-LAL,ANALY_ARC,ANALY_GLASGOW"
    excludeSite += ",ANALY_TRIUMF_CVMFS,ANALY_CERN_XROOTD,ANALY_CAM,ANALY_GRIF-IRFU,ANALY_DRESDEN"
    excludeSite += ",ANALY_TR-10-ULAKBIM,ANALY_GR-01-AUTH,ANALY_LPC"
    excludeSite += ",ANALY_ZA-UJ,ANALY_SFU_bugaboo"
    submit   = "prun "# --excludedSite="+excludeSite
    submit  += " --useSiteGroup=1"   # submit only to Tier1 + alpha (most reliable) sites
    submit  += " --maxCpuCount=43200"  # prevent "Looping job killed by pilot" after 2 hours w/o printout - increase to 8 hours
#   submit  += " --memory=4000"   # to exclude sites that don't allow >4000MB jobs (default mostly 3000; maybe required with 64-bit and large workspace)
    submit  += " --useChirpServer=voatlas92.cern.ch"
    submit += " --rootVer=6.04/16"
    submit += " --cmtConfig=x86_64-slc6-gcc49-opt"
#   submit  += " --destSE=UKI-SOUTHGRID-RALPP_LOCALGROUPDISK"
    if not "-u" in opt:
        vomsGroup= "phys-higgs"
        submit += " --official --voms=atlas:/atlas/"+vomsGroup+"/Role=production"
    job      = "./"+jobname+"_job_grid"
    jobsplit = 100   # Alastair Dewhurst says 300-500 jobs will fill a smaller site (17/5/2011)
    sleepTime=  30   # Give the Panda monitor enough time to see the site is full with out last submission (Alastair Dewhurst says it updates every 2mins)

    if "-u" in opt:
        user    = pwd.getpwuid(os.geteuid())[0]    # os.getlogin() doesn't always work (eg. with nohup)
        dsprefix= "user."+user
    else:
        dsprefix= "group."+vomsGroup+".comb"

else:
    job      = "./"+jobname+"_job"
    if   "-r" in opt:
        submit = ""                           # run directly
    elif "-s" in opt:
        submit = "~adye/bin/job -O"           # use my background job command
    else:
        submit = "bsub -q 1nh"                # CERN lxplus
        job    = jobname+"_job"
    sleepTime= 0
