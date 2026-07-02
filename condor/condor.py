#!/usr/bin/env python
##
## @package condor.py
## @file condor.py
## @brief python methods for running condor
##
## @brief Created by:
## @author John R. Peterson (Purdue)
##
## @brief Modified by:
## @author Emily Grace (Purdue)
## @author Nathan Todd (Purdue)
## @author En-Hsin Peng (Purdue)
## @author Glenn Sembroski (Purdue)
## @author Jim Chiang (SLAC)
## @author Jeff Gardner (Google)
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.

import os
import subprocess

## Condor method to setup directories
def initEnvironment(self):
    self.dagfile=open(self.workDir+'/dag_'+self.observationID+'.dag','w')
    if not os.path.exists(self.workDir+'/logs'):
        os.makedirs(self.workDir+'/logs')
    if not os.path.exists(self.workDir+'/errors'):
        os.makedirs(self.workDir+'/errors')
    if not os.path.exists(self.workDir+'/output'):
        os.makedirs(self.workDir+'/output')

## Condor method to create the submit jobs (trim,raytrace,e2adc)
def writeSubmit(self, job, jobName, fid='none', ckpt=0):
    if self.grid != 'condor':
        return
    assert job in ['raytrace', 'trim', 'e2adc']
    assert 'universe' in self.grid_opts
    universe = self.grid_opts['universe']
    submitfile=open(jobName+'.submit','w')
    submitfile.write('executable = %s/%s\n' % (self.binDir,job))
    submitfile.write('initialdir = %s\n' % self.workDir)
    submitfile.write('Universe = %s\n' % universe)
    submitfile.write('Input = %s.pars\n' % jobName)

    if job == 'raytrace' or job == 'e2adc':
        submitfile.write('stream_output = True\n')
        submitfile.write('Error = errors/error_%s.error\n' % fid)
    else:
        submitfile.write('Log = logs/log_%s.log\n' %jobName)
        submitfile.write('Output = output/out_%s.out\n' %jobName)
        submitfile.write('Error = errors/error_%s.error\n' %jobName)

    if universe == 'vanilla':
        submitfile.write('periodic_release=true\n')
        submitfile.write('should_transfer_files = YES\n')
    else:
        submitfile.write('should_transfer_files = IF_NEEDED\n')
    submitfile.write('when_to_transfer_output = ON_EXIT\n')
    submitfile.write('notification = NEVER\n')

    instrument=self.instrDir.split("/")[-1]
    fidfilt=fid.replace(self.observationID,self.observationID+'_f'+self.filt)
    if job == 'raytrace':
        submitfile.write('Log = logs/log_r_%s.log\n' % fid)
        submitfile.write('Output = output/out_r_%s.out\n' % fid)
        submitfile.write('transfer_input_files = \\\n')
        submitfile.write('airglowscreen_%s.fits.gz, \\\n' % self.observationID)
        submitfile.write('fea_%s_0.txt, \\\n' % self.observationID)
        submitfile.write('fea_%s_1.txt, \\\n' % self.observationID)
        submitfile.write('fea_%s_2.txt, \\\n' % self.observationID)
        #Pick up the number ot Atmosphere layers from the
        #  'atmosphere_"self.observationID".pars' file which was made by the 
        #atmosphere program.  The default is 8, set in atmosphere.main but I think 
        #it can be overridden by a 'natmospherefile' comand in the physics command 
        #file.
        atmFile=('atmosphere_'+self.observationID+'.pars')
        afile=open(atmFile,'r')
        for line in afile:
            lstr=line.split()
            if "natmospherefile " in line:
                natmospherefile=lstr[1]
            if "cloudfile " in line:
                submitfile.write('cloudscreen_%s_%s.fits.gz, \\\n' % (self.observationID,lstr[1]))
	       
        for layer in range(int(natmospherefile)):
            for f in ['coarsep','coarsex','coarsey','fineh','finep',
                      'largep','largex','largey','mediumh','mediump','mediumx',
                       'mediumy']:
                submitfile.write('atmospherescreen_%s_%d_%s.fits.gz, \\\n' % 
                                                  (self.observationID,layer,f))

        # At this point we want to transfer over all of the insturment 
        # directory files. This enables different numbers of fiters and 
        # optics to be all copied.
        instDirList=os.listdir(self.instrDir)
        for file in instDirList :
            absoluteFilePath = self.instrDir + '/' + file
            submitfile.write('%s,\\\n' %  absoluteFilePath)

        # And do the same for all the files in data/sky. This includes 
        # sed_flat.txt and files needed for background.
        skyDirList=os.listdir(self.dataDir+'/sky')
        for file in skyDirList :
            absoluteFilePath = self.dataDir + '/sky/' + file
            submitfile.write('%s,\\\n' %  absoluteFilePath)

        for f in ['climate','earth','exosphere_1000','exosphere_1200',
                  'exosphere_1400','exosphere_1600','exosphere_500',
                  'exosphere_600','exosphere_700','exosphere_800',
                  'exosphere_900','h2ocs','midlatitude_summer',
                  'midlatitude_winter','o2cs','o3cs','population',
                  'subarctic_summer','subarctic_winter','tropical',
                  'us_standard','us_standard_trace_1','us_standard_trace_2',
                  'us_standard_trace_3','us_standard_trace_4','us_standard_trace_5',
                  'wind']:
            submitfile.write('%s/atmosphere/%s.txt, \\\n' % (self.dataDir,f))

        materialDirList=os.listdir(self.dataDir+'/material')
        for file in materialDirList:
            absoluteFilePath = self.dataDir + '/material/' + file
            submitfile.write('%s,\\\n' %  absoluteFilePath)

        submitfile.write('%s/version, \\\n'% (self.binDir))

        # ################################################################
        # This now gets a little tricky. We don't actually finish the 
        # transfer_input_file list here. We still need to add whatever SEDs 
        # files this raytrace job will need. However we won't know that until
        # the trim programs make up the "trimmed" catalogs. We will also need 
        # to add the "Queue 1" command to the end of the submit file (needs to
        # be the last line). All that is done in the  condor/chip::posttrim() 
        # function that is called by condor after the trim program completes but 
        # before the raytrace submissions are made. 
        # ################################################################
        submitfile.write('tracking_%s.pars ' % self.observationID)

        if ckpt>0:
            submitfile.write(', %s_e_%s_ckptdt_%d.fits.gz' % 
                                                   (instrument,fidfilt,ckpt-1))
            submitfile.write(', %s_e_%s_ckptfp_%d.fits.gz ' % 
                                                   (instrument,fidfilt,ckpt-1))
        # the chip::posttrim() function will add the SED files and the "Queue 1"
        # command. See comment above
        
    else:
        if job == 'trim':
            submitfile.write('transfer_input_files = \\\n')
            instDirList=os.listdir(self.instrDir)
            for file in instDirList :
                absoluteFilePath = self.instrDir + '/' + file
                submitfile.write('%s,\\\n' %  absoluteFilePath)
#            submitfile.write('transfer_input_files =  %s/focalplanelayout.txt' % self.instrDir)
            for line in open('catlist_'+self.observationID+'.pars'):
                submitfile.write(', %s' % line.split()[2])
                submitfile.write('\n')
        elif job == 'e2adc':
            submitfile.write('Log = logs/log_e_%s.log\n' % fid)
            submitfile.write('Output = output/out_e_%s.out\n' % fid)
            submitfile.write('transfer_input_files = %s/segmentation.txt, ' % self.instrDir)
            submitfile.write('%s/focalplanelayout.txt, ' % self.instrDir)
            submitfile.write('%s_e_%s.fits.gz\n' % (instrument,fidfilt))
        submitfile.write('Queue 1\n')
    submitfile.close()

def writeTrimDag(self,jobName,tc,nexp):
    checkpoint=self.grid_opts.get('checkpoint', 12)
    self.dagfile.write('JOB %s %s/%s.submit\n' % (jobName,self.workDir,jobName))
    #condor/chip.postTrim needs sedDir.
    self.dagfile.write('SCRIPT POST %s %s/condor/chip posttrim %s %d %d %s %s %s %d\n' 
                       % (jobName,self.phosimDir,self.observationID,tc,
                          nexp,self.dataDir,self.workDir,self.sedDir, checkpoint))
    self.dagfile.write('RETRY %s 20\n' % (jobName))
    writeSubmit(self,'trim',jobName)

def writeRaytraceDag(self,cid,eid,tc,run_e2adc):
    checkpoint=self.grid_opts.get('checkpoint', 0)   #No checkpointing is default.
    observationID=self.observationID
    fid=observationID + '_' + cid + '_' + eid
    fidfilt=observationID + '_f' + self.filt + '_' + cid + '_' + eid
    instrument=self.instrDir.split("/")[-1]
    if os.path.exists(self.outputDir+'/'+instrument+'_'+fidfilt+'.tar'):
        os.remove('raytrace_'+fid+'.pars')
        return

    if not os.path.exists(instrument+'_e_'+fidfilt+'.fits.gz'):

        #This is a little tricky. If we want no checkpointing (checkpoint = 0),
        #then just do 1 raytrace_8_0.para (.submit)

        firstckpt=0
        for ckpt in range(checkpoint+1):
            if os.path.exists(instrument+'_e_'+fidfilt+'_ckptdt_'+str(ckpt)+'.fits.gz') and os.path.exists(instrument+'_e_'+fidfilt+'_ckptfp_'+str(ckpt)+'.fits.gz'):
                firstckpt=ckpt+1

        # if no checkpinitng: firstckpt and checkpoint is 0 and we only make the first
        # raytrace*_0.pars,  and .submit files and add their commands to the dagfile 
        for ckpt in range(firstckpt,checkpoint+1):
            fidckpt=fid+'_'+str(ckpt)
            self.dagfile.write('JOB raytrace_%s %s/raytrace_%s.submit\n' % (fidckpt,self.workDir,fidckpt))
            self.dagfile.write('RETRY raytrace_%s 20\n' % (fidckpt))
            if ckpt==firstckpt:
                self.dagfile.write('PARENT trim_%s_%d CHILD raytrace_%s\n' % (observationID,tc,fidckpt))
            else:
                self.dagfile.write('PARENT raytrace_%s_%d CHILD raytrace_%s\n' % (fid,ckpt-1,fidckpt))

            if ckpt==checkpoint:
                if run_e2adc:
                    self.dagfile.write('SCRIPT POST raytrace_%s %s/condor/chip lastraytracecleanup %s %s %s %s %d %s %s\n' %
                            (fidckpt,self.phosimDir,observationID,self.filt,cid,eid,ckpt,self.workDir, instrument))
                else:
                    self.dagfile.write('SCRIPT POST raytrace_%s_%d %s/condor/chip postraytrace %s %s %s %s %d %s %s %s\n' %
                            (fid,checkpoint,self.phosimDir,observationID,self.filt,cid,eid,checkpoint,self.workDir,self.outputDir,instrument))
            else:
                self.dagfile.write('SCRIPT POST raytrace_%s %s/condor/chip raytracecleanup %s %s %s %s %d %s %s %s\n' %
                        (fidckpt,self.phosimDir,observationID,self.filt,cid,eid,ckpt,self.workDir, self.outputDir, instrument))

            writeSubmit(self,'raytrace','raytrace_'+fidckpt,fid,ckpt)
            pfile=open('raytrace_'+fidckpt+'.pars','w')
            pfile.write(open('raytrace_'+fid+'.pars').read())
            pfile.write("checkpointcount %d\n" % ckpt)
            pfile.write("checkpointtotal %d\n" % checkpoint)
            pfile.close()
    if run_e2adc:
        self.dagfile.write('JOB e2adc_%s %s/e2adc_%s.submit\n' % (fid,self.workDir,fid))
        self.dagfile.write('RETRY e2adc_%s 20\n' % fid)
        self.dagfile.write('SCRIPT POST e2adc_%s %s/condor/chip poste2adc %s %s %s %s %s %s %s %s\n' %
                      (fid,self.phosimDir,observationID,self.filt,cid,eid,self.outputDir,self.instrDir,self.workDir,instrument))
        if not os.path.exists(instrument+'_e_'+fidfilt+'.fits.gz'):
            self.dagfile.write('PARENT raytrace_%s_%d CHILD e2adc_%s\n' % (fid,checkpoint,fid))
        else:
            self.dagfile.write('PARENT trim_%s_%d CHILD e2adc_%s\n' % (observationID,tc,fid))
        writeSubmit(self,'e2adc','e2adc_'+fid,fid)
    else:
        if os.path.exists(instrument+'_e_'+fidfilt+'.fits.gz'):
            command=('%s/condor/chip postraytrace %s %s %s %s %d %s %s %s' %
                    (self.phosimDir,observationID,self.filt,cid,eid,checkpoint,self.workDir,self.outputDir,instrument))
            if subprocess.call(command, shell=True) != 0:
                raise RuntimeError("Error running %s" % command)
    os.remove('raytrace_'+fid+'.pars')

def submitDag(self):
    self.dagfile.close()
    command='condor_submit_dag -MaxPre 8 -MaxPost 8 -AutoRescue 0 dag_'+self.observationID+'.dag'
    if subprocess.call(command, shell=True) != 0:
        raise RuntimeError("Error running %s" % command)

