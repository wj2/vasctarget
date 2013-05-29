import VisionEgg, pygame, datetime, time, math, StringIO, sys
from VisionEgg.Core import *; from VisionEgg.MoreStimuli import *;
from Dimstim.Core import *

def sparsenoise(preexpSec, postexpSec, orioff,
                ncellswide, ncellshigh, regionwidthDeg, regionheightDeg,
                xi, yi, ori, widthDeg, heightDeg, brightness, colour, antialiase, bgbrightness, sweeptimeMsec, postsweepMsec,
                varlist, nruns, shuffleRuns, blankSweep, shuffleBlankSweeps):

	STIMTYPE = 'Sparse noise'

	#######################################################################################
	# Define parameter lists
	paramlist = ['preexpSec','postexpSec','orioff','ncellswide','ncellshigh','regionwidthDeg','regionheightDeg',
                 'xi','yi','ori','widthDeg','heightDeg','brightness','colour','antialiase','bgbrightness','sweeptimeMsec','postsweepMsec',
                 'varlist','nruns','shuffleRuns','blankSweep','shuffleBlankSweeps'] # list of all parameters passed to this function
	NVSparamlist = paramlist[:] # copy it
	map(NVSparamlist.remove,['preexpSec','postexpSec', 'colour', 
	                         'varlist','shuffleRuns']) # list of all user specifiable experiment paramaters which have a field in NVS header (subset of paramlist)
	dynvarnames = ['xi','yi','ori','widthDeg','heightDeg','brightness','colour','antialiase','bgbrightness','sweeptimeMsec','postsweepMsec'] # list of all potentially dynamic variables (subset of paramlist)

	# Do some initial error checking
	for var in dynvarnames:
		if type(eval(var)) in [list,tuple] and var not in varlist and var <> 'colour' :
			raise RuntimeError('Variable \'%s\' has been assigned multiple values, but it hasn\'t been entered into the varlist' %var)
	for var in varlist:
		if type(eval(var)) not in [list,tuple]:
			raise RuntimeError('Variable \'%s\' has been entered into the varlist, yet it isn\'t a list. Add [brackets] around it or remove it from varlist' %var)

	#######################################################################################
	# Load attribs from Dimstim config file
	origDeg = eval(DC.get('Manbar1','POS')) # deg, wrt screen center
	if orioff == None: # if it doesn't exist, get value from Manbar1 setting
		orioff = eval(DC.get('Manbar1','ORI')) # deg
	if widthDeg == None: # if it doesn't exist, get value from Manbar1 setting
		widthDeg = eval(DC.get('Manbar1','SIZE'))[0] # deg
	elif widthDeg == [None]: # if it doesn't exist, get value from Manbar1 setting, leave as list
		widthDeg = [ eval(DC.get('Manbar1','SIZE'))[0] ] # deg
	if heightDeg == None: # if it doesn't exist, get value from Manbar1 setting
		heightDeg = eval(DC.get('Manbar1','SIZE'))[1] # deg
	elif heightDeg == [None]: # if it doesn't exist, get value from Manbar1 setting, leave as list
		heightDeg = [ eval(DC.get('Manbar1','SIZE'))[1] ] # deg
	noWidthOverlap = 0
	if regionwidthDeg == None: # if it doesn't exist, set to be multiple of widthDeg (no overlap)
		if type(widthDeg) is not list:
			regionwidthDeg = ncellswide*widthDeg #deg
		elif type(widthDeg) is list:
			regionwidthDeg = ncellswide*widthDeg[0] #deg
		noWidthOverlap = 1
	noHeightOverlap = 0
	if regionheightDeg == None: # if it doesn't exist, set to be multiple of heightDeg (no overlap)
		if type(heightDeg) is not list:
			regionheightDeg = ncellshigh*heightDeg #deg
		elif type(heightDeg) is list:
			regionheightDeg = ncellshigh*heightDeg[0] #deg
		noHeightOverlap = 1

	# Generate sweep table and sweeplist
	varvals={} # init a dictionary that will contain variable values
	for var in varlist:
		varvals[var]=eval(var) # generate a dictionary with var:val entries, to pass to buildSweepTable
	(sweepTable,dimlist,sweeplist,sweeptabletext) = buildSweepTable(varlist,varvals,nruns,shuffleRuns,blankSweep,shuffleBlankSweeps,makeSweepTableText=0) # passing varlist by reference, dim indices end up being modified
        #if EYELINKINSTALLED:
        #    EYELINK.sendMessage("sweeplist ="+str(sweeplist))
        #print "X"*10, "sweeplist ="+str(sweeplist)
	nsweeps = len(sweeplist)

	# Do time and space conversions of applicable static parameters
	orig=[[],[]] # XY origin, in pix from bottom left, init as a 2 item list
	orig[0] = SCREENWIDTH/2.0 + deg2pix(origDeg[0]) # x pix
	orig[1] = SCREENHEIGHT/2.0 + deg2pix(origDeg[1]) # y pix
	preexp  = sec2intvsync(preexpSec); preexpSec = vsync2sec(preexp) # pre-experiment duration, in vsyncs; update requested value
	postexp = sec2intvsync(postexpSec); postexpSec = vsync2sec(postexp) # post-experiment duration, in vsyncs; update requested value
	blanksweeptime = sec2intvsync(blankSweep[1]); blankSweep = (blankSweep[0], vsync2sec(blanksweeptime)) # blank sweep duration, in vsyncs, update requested value in blankSweep tuple
	regionwidth  = deg2pix(regionwidthDeg) # grid width, in pix
	regionheight = deg2pix(regionheightDeg) # grid height, in pix
	# dec or inc (whichever's closer) to an integer number of screen pixels wide and high per grid element
	regionwidth  = quantizeSpace(regionwidth,ncellswide) # returns an int
	regionwidthDeg = pix2deg(regionwidth) # update requested value
	regionheight = quantizeSpace(regionheight,ncellshigh) # returns an int
	regionheightDeg = pix2deg(regionheight) # update requested value
	gwspacing = regionwidth/float(ncellswide) # this should divide evenly now
	ghspacing = regionheight/float(ncellshigh)
	#print 'ncellswide is',ncellswide
	#print 'ncellshigh is',ncellshigh
	#print 'regionwidth is',regionwidth
	#print 'regionheight is',regionheight
	#print 'gwspacing is',gwspacing
	#print 'ghspacing is',ghspacing

	# Do time and space conversions of applicable potentially dynamic variables, convert from sweepTable to sweeptable
	sweeptable = sweepTable.copy() # init converted sweeptable to be a copy of unconverted sweepTable
	for var in dynvarnames: # for all potentially dynamic variables (ie, all those that could potentially be in sweeptable)
		if var in ['widthDeg','heightDeg']:
			varDeg = var
			var=var.replace('Deg','') # remove the units from the var name
			if varDeg in sweeptable: # if var is in sweeptable
				temp=[]
				for deg in sweeptable[varDeg]:
					temp.append(deg2pix(deg)) # convert to pix, can't quantize cuz it's changing on sweeps
				sweeptable[var] = temp # add new entry to sweeptable with these units
				del sweeptable[varDeg] # remove varDeg entry from sweeptable
			else: # varDeg isn't in sweeptable
				exec(var+'=deg2pix('+varDeg+')') # convert and init it
				if var == 'width' and noWidthOverlap == 1:
					width = gwspacing # quantize it to the grid spacing
					print 'width is', width
					widthDeg = pix2deg(width) # update requested value
				elif var == 'height' and noHeightOverlap == 1:
					height = ghspacing # quantize it to the grid spacing
					print 'height is', height
					heightDeg = pix2deg(height) # update requested value
		elif var in ['sweeptimeMsec','postsweepMsec']:
			varMsec = var
			var=var.replace('Msec','') # remove the units from the var name
			if varMsec in sweeptable: # if var is in sweeptable
				vsynctemp=[]; msectemp = []
				for msec in sweeptable[varMsec]:
					vsynctemp.append(msec2intvsync(msec)) # convert to integer vsyncs
					msectemp.append(vsync2msec(vsynctemp[-1])) # convert back to msec to update original
				sweeptable[var] = vsynctemp # add new entry to sweeptable with these units
				sweepTable[varMsec] = msectemp # update varMsec entry in sweepTable
				del sweeptable[varMsec] # remove varMsec entry from sweeptable
			else: # varMsec isn't in sweeptable
				exec(var+'=msec2intvsync('+varMsec+')') # convert to integer vsyncs and init it
				exec(varMsec+'=vsync2msec('+var+')') # update requested value
	'''
	# Print units converted sweeptable
	print 'units converted sweeptable:'
	for var in sweeptable:
		print var+' =', sweeptable[var]
	print
	'''
	# Calculate grids
	if type(ori) is not list: # get ori as a list first
		orilist = [ori]
	else:
		orilist = ori
	grid={} # make it a dictionary, indexed into using ori, xi and yi
	for tempori in orilist:
		tempgrid=[]
		for i in range(0,ncellswide):
			tempgrid.append([])
			for j in range(0,ncellshigh):
				tempgrid[i].append([])
				x0 = (0.5+i)*gwspacing-regionwidth/2.0
				y0 = (0.5+j)*ghspacing-regionheight/2.0
				r = math.sqrt(x0**2+y0**2)
				if x0 == 0.0:
					if y0 > 0.0: # on the +ve y-axis
						theta = math.pi/2.0
					elif y0 < 0.0: # on the -ve y-axis
						theta = 3/2.0*math.pi
					elif y0 == 0.0: # at the origin, doesn't matter what theta is, set it to zero anyway
						theta = 0.0
				else: # x0 != 0, so can safely do arctan(y0/x0)
					theta = abs(math.atan(y0/float(x0))) # in rad
					if x0 > 0.0:
						if y0 > 0.0: # in the 1st quadrant
							pass
						elif y0 < 0.0: # in the 4th quadrant
							theta = 2*math.pi - theta
						elif y0 == 0.0: # on +ve x-axis
							theta = 0.0
					elif x0 < 0.0:
						if y0 > 0.0: # in the 2nd quadrant
							theta = math.pi - theta
						elif y0 < 0.0: # in the 3rd quadrant
							theta = math.pi + theta
						elif y0 == 0.0: # on -ve x-axis
							theta = math.pi
				x = r*math.cos(theta+math.radians(tempori+orioff))+orig[0] # grid oriented with offset and centered
				y = r*math.sin(theta+math.radians(tempori+orioff))+orig[1] # grid oriented with offset and centered
				tempgrid[i][j] = (x,y)
		grid[tempori]=tempgrid # set tempgrid as the dictionary entry for tempori

        print grid
	# Calculate exptimeSec
	if 'sweeptime' in sweeptable: # if sweeptime varies over sweeps, sum up over all sweeps in sweeplist
		nsweepvsyncs = sum( [ sweeptable['sweeptime'][sweepi] for sweepi in sweeplist if sweepi is not None ] ) # all but the blank sweeps
	else: # sweeptime is constant over all sweeps excluding blank sweeps
		nsweepvsyncs = sum( [ sweeptime for sweepi in sweeplist if sweepi is not None ] ) # all but the blank sweeps
	nsweepvsyncs += sum( [ blanksweeptime for sweepi in sweeplist if sweepi is None ] ) # add blank sweeps, displayed for blanksweeptime duration
	if 'postsweep' in sweeptable: # if postsweep varies over sweeps, sum up over all sweeps in sweeplist
		npostsweepvsyncs = sum( [ sweeptable['postsweep'][sweepi] for sweepi in sweeplist if sweepi is not None ] ) # all but the blank sweeps, blank sweeps have no postsweep
	else: # postsweep is constant over all sweeps
		npostsweepvsyncs = sum( [ postsweep for sweepi in sweeplist if sweepi is not None ] ) # all but the blank sweeps, blank sweeps have no postsweep
	exptimeSec = vsync2sec(nsweepvsyncs + npostsweepvsyncs) # in sec
	# there's no delay to do pre-sweep calculations for sparse noise
	print 'Experiment time will be',nicetime(exptimeSec,6)
	print

	# Set NVS header and text header
	hdr = Header(exptimeSec)
	# set internal (not user-specified in the experiment) static parameters that have a field in NVS header
	hdr.NVS[STT_ORIGX]        = origDeg[0] # X origin wrt screen center (deg)
	hdr.NVS[STT_ORIGY]        = origDeg[1] # Y origin wrt screen center (deg)
	hdr.NVS[STT_STS]          = 11 # stimulus code for sparse noise
	hdr.NVS[STT_EYE]          = EYE
	hdr.NVS[STT_TOTAL_SWEEPS] = nsweeps # num sweeps in entire experiment
	# set everything else (user-specified static parameters) that's applicable in NVS header
	NVSparamdict={} # init a dictionary that will contain variable values
	for var in NVSparamlist:
		NVSparamdict[var]=eval(var) # generate a dictionary with var:val entries, to pass to buildNVS
 	buildNVS(hdr,NVSparamdict,varlist,dimlist)
	# set text header
	hdr.text = 'STIMTYPE = \''+STIMTYPE+'\'\n'+ \
	           'EXPERIMENT = \''+EXPERIMENT+'\'\n'+ \
	           'SCREENWIDTHCM = '+str(SCREENWIDTHCM)+'\n'+ \
	           'SCREENHEIGHTCM = '+str(SCREENHEIGHTCM)+'\n'+ \
	           'SCREENDISTANCECM = '+str(SCREENDISTANCECM)+'\n'+ \
	           'SCREENWIDTH = '+str(SCREENWIDTH)+'\n'+ \
	           'SCREENHEIGHT = '+str(SCREENHEIGHT)+'\n'+ \
	           'REFRESHRATE = '+str(REFRESHRATE)+'\n'+ \
	           'PIXPERCM = '+str(PIXPERCM)+'\n'+ \
	           'GAMMA = '+str(GAMMA)+'\n'+ \
	           'GAMMAOFFSET = '+str(GAMMAOFFSET)+'\n'+ \
	           'EYE = '+str(EYE)+'\n'+ \
	           'FIXATION = '+str(FIXSPOT)+'\n'+ \
	           'origDeg = '+str(origDeg)+'\n'
	for param in paramlist:
		if type(eval(param)) is types.StringType:
			hdr.text += param+' = \''+eval(param)+'\'\n' # add quotes around the value
		else:
			hdr.text += param+' = '+str(eval(param))+'\n'
	hdr.text += 'dimlist = '+str(dimlist)+'\n' # not absolutely necessary, but handy, and not very long, so why not...
	hdr.text += 'nsweeps = '+str(nsweeps)+'\n'
	longsweepliststring = 'sweeplist = '+str(sweeplist).replace(' ','')+'\n' # remove spaces to save space
	if len( hdr.text + longsweepliststring ) <= TEXTLENGTH: # if adding the long sweeplist string to text header won't surpass the text header limit
		hdr.text += longsweepliststring # then add the long sweeplist string

	# Print the NVS and text headers
	#hdr.printNVS(); print; # print the NVS header
	print 'Text header length is', len(hdr.text)
	hdr.printText(); # print the text header

	# Init OpenGL graphics screen
	screen = get_default_screen()
	sp = screen.parameters

	# Create an instance of the Target2D class
	target = Target2D(anchor = 'center',
	                  on     = 0) # keep it off until first sweep starts
	tp = target.parameters
        if 'bgbrightness' in sweeptable:
            background = Target2D(anchor = 'center',
                                  size = (SCREENWIDTH,SCREENHEIGHT),
                                  position = (SCREENWIDTH/2.0,SCREENHEIGHT/2.0),
                                  color = (bgbrightness[0],bgbrightness[0],bgbrightness[0],1.0))

        else:
            background = Target2D(anchor = 'center',
                                  size = (SCREENWIDTH,SCREENHEIGHT),
                                  position = (SCREENWIDTH/2.0,SCREENHEIGHT/2.0),
                                  color = (bgbrightness,bgbrightness,bgbrightness,1.0))
	bg = background.parameters

	# Create a Viewport instance
	#viewport = Viewport(screen=screen, stimuli=[background,fspot,target])
	viewport = Viewport(screen=screen, stimuli=[background,target])

	# Draw viewport once in advance, fixes weird lag on first sweep
	tp.on = 1
	viewport.draw()
	tp.on = 0

	# Shows fixation spot for drift correction, if EYELINKINSTALLED
	screen.clear()
        viewport.draw()
	swap_buffers() # returns on next vsync pulse from video card

        # If Eyelink installed, run a pre-stimulus drift correction, start recording...
        if EYELINKINSTALLED:
            #error = EYELINK.doDriftCorrect(int(SCREENWIDTH/2.0),int(SCREENHEIGHT/2.0),0,0)
            #if (error != 27 or error != -1): # no problem with drift correction
                error = EYELINK.startRecording(1,1,0,0) # record samples and events to EDF, but don't send over the link
                if error: raise RuntimeError('Error starting recording') # quit
                msecDelay(50) # ensure recording started before starting stimulus
                EYELINK.sendMessage('Stimulus started')
            #else:
            #    raise RuntimeError('Eyelink drift correction unsatisfactory');

	# Init DT board and send the entire stimulus header, including NVS and text headers, plus extra info before and after that SURF needs
	exec('DT.InitBoard()')*DTBOARDINSTALLED
	#exec('headerchksum = hdr.broadcast()')*DTBOARDINSTALLED
        # set this immediately so that we know when condition zero occurs
        exec('DT.postInt16NoDelay(BLANKSWEEP)')*DTBOARDINSTALLED # sweep DOUT

	# Init vars for main loop
	nvsyncsdisplayed = 0 # 1 based, init as 0
	vsynctimer = VsyncTimer()
	for var in sweeptable: # init local copies of vars in sweeptable to their vals at first sweep in sweeplist
		exec(var + '= sweeptable[\'' + var + '\'][sweeplist[0]]') # ie set var = sweeptable[var][sweeplist[0]], e.g. ori=sweeptable['ori'][sweeplist[0]]

	bg.color = (bgbrightness,bgbrightness,bgbrightness,1.0) # set bg colour now that bgbrightness is init'd, do this now so it's correct for the pre-exp delay

	# Do pre-experiment delay
	(blah,quit,pause,pausehappened,resumehappened)=staticScreen(screen=screen,
                                                                    viewport=viewport,
                                                                    nvsyncs=preexp)

        ## If Eyelink installed, run a pre-stimulus drift correction, start recording...
        #if EYELINKINSTALLED:
        #    #error = EYELINK.doDriftCorrect(int(SCREENWIDTH/2.0),int(SCREENHEIGHT/2.0),0,0)
        #    #if (error != 27 or error != -1): # no problem with drift correction
        #        error = EYELINK.startRecording(1,1,1,0) # record samples and events to EDF, AND send over the link
        #        if error: raise RuntimeError('Error starting recording') # quit
        #        msecDelay(50) # ensure recording started before starting stimulus
        #        EYELINK.sendMessage('Stimulus started')
        #    #else:
        #    #    raise RuntimeError('Eyelink drift correction unsatisfactory');

	expBeganDateTime = datetime.datetime.now()
	expBeganTime = time.clock() # precision timestamp

	# Main loop
        #print "M"*100
        #print sweeptable
	for (sweeplisti,sweepi) in enumerate(sweeplist):

		if sweepi is None: # do a blank sweep
			tp.on = 0 # turn off the stimulus, leave all other parameters unchanged
			postval = BLANKSWEEP # posted to DT port to indicate a blank sweep
			nvsyncs = blanksweeptime # this many vsyncs for this sweep
			npostvsyncs = 0 # this many post-sweep vsyncs for this sweep, blank sweeps have no post-sweep delay
		else: # not a blank sweep
			tp.on = 1 # ensure stimulus is on
			postval = sweepi # sweep index number will be posted to DT port
			for var in sweeptable: # load dynamic variables for this sweepi from sweeptable
				exec(var + '= sweeptable[\'' + var + '\'][sweepi]') # ie set var = sweeptable[var][sweepi], e.g. ori=sweeptable['ori'][sweepi]
			nvsyncs = sweeptime # this many vsyncs for this sweep
			npostvsyncs = postsweep # this many post-sweep vsyncs for this sweep

		# Update target parameters
		tp.position      = grid[ori][xi][yi]
		tp.orientation   = orioff + ori
		tp.size          = (width,height)
		tp.color         = (colour[0]*brightness,colour[1]*brightness,colour[2]*brightness,1.0)
		tp.anti_aliasing = antialiase

		# Update background parameters
                white = [1,1,1]
                #print sweeptable
		bg.color = (bgbrightness,bgbrightness,bgbrightness,1.0)

		# Set sweep bit high, do the sweep
		#exec('DT.setBitsNoDelay(SWEEP)')*DTBOARDINSTALLED # set sweep bit high, no delay
                exec('DT.postInt16NoDelay(postval+SWEEP)')*DTBOARDINSTALLED # post value to port, no delay
		for vsynci in xrange(0,nvsyncs): # range depends on if this is a blank sweep or not
			for eventi in pygame.event.get(): # for all events in the event queue
				if eventi.type == pygame.locals.KEYDOWN:
					if eventi.key == pygame.locals.K_ESCAPE:
						quit = 1
					if eventi.key == pygame.locals.K_PAUSE:
						pause = int(not pause) # toggle pause
						#pausehappened = 1
			# check for broken fixation
			#fbit = 1 # if no DT340 board, will not pause by default
			#exec('fbit=DT.getByte()')*DTBOARDINSTALLED				
			#if not(fbit & FIXN): # fixation broken, pause stimulus
                        #        pause = 1
			#if quit:
			#	break # out of vsync loop
                        if EYELINKINSTALLED:
                                EYELINK.sendMessage("%d" %(postval)+chr(9)+"%d"%(vsynci))
			screen.clear()
			viewport.draw()
			swap_buffers() # returns on next vsync pulse from video card
			vsynctimer.tick()
			nvsyncsdisplayed += 1 # increment

		# Sweep's done, turn off the target, clear sweep bit low, do the postsweep delay
		tp.on = 0
		#exec('DT.clearBitsNoDelay(SWEEP)')*DTBOARDINSTALLED # clear sweep bit low, no delay
		(vsynctimer,quit,pause,pausehappened,resumehappened)=staticScreen(screen=screen,
		                                                     viewport=viewport,
		                                                     nvsyncs=npostvsyncs, # depends on if this is a blank sweep or not
		                                                     vsynctimer=vsynctimer,
		                                                     quit=quit,
		                                                     pause=pause,
		                                                     pausehappened=pausehappened,
                                                                     resumehappened=resumehappened)
		if quit:
			sweeplisti -= 1 # dec for accurate count of how many sweeps were displayed fully
			break # out of sweep loop

	expEndedTime = time.clock() # precision timestamp
	expEndedDateTime = datetime.datetime.now()
	exec('DT.postInt16NoDelay(0)')*DTBOARDINSTALLED # clear the value, no delay

	# Do post-experiment delay, delay for postexp time less post-sweep time of last sweep
	(blah,quit,pause,pausehappened,resumehappened)=staticScreen(screen=screen,
                                                                    viewport=viewport,
                                                                    nvsyncs=postexp-npostvsyncs,
                	                                            quit=quit,
                	                                            pause=pause,
                                                                    pausehappened=pausehappened,
                                                                    resumehappened=resumehappened)

	# Send checksum
	#exec('sweepchksum = DT.getChecksum()')*DTBOARDINSTALLED
	#exec('DT.postInt16NoDelay(sweepchksum)')*DTBOARDINSTALLED # post the value, no delay
	#exec('DT.toggleBits(DATA)')*DTBOARDINSTALLED # toggle data bit, delay (triggers SURF)

	# Dimstim is no longer running, close DT board
	exec('DT.clearBitsNoDelay(0x00ffffff)')*DTBOARDINSTALLED # clear all bits (including run bit) low, no delay
	exec('DT.CloseBoard()')*DTBOARDINSTALLED

	# Print messages to VisionEgg log and to screen
	buf = StringIO.StringIO()
	exec("print >> buf, 'Sweep table:'")*(sweeptabletext != None)
	exec("print >> buf, sweeptabletext")*(sweeptabletext != None)
	print >> buf, 'Stimulus type: _______________', STIMTYPE
	print >> buf, 'Experiment: __________________', EXPERIMENT.replace('\\\\','\\')
	print >> buf, 'Experiment began: ____________', expBeganDateTime
	print >> buf, 'Experiment ended: ____________', expEndedDateTime
	print >> buf, 'Predicted experiment time: ___ %s' %(nicetime(exptimeSec,6))
	print >> buf, 'Actual experiment time: ______ %s' %(nicetime(expEndedTime - expBeganTime,6))
	print >> buf, 'Number of runs: ______________ %s' %nruns
	print >> buf, 'Blank sweep every: ___________ %s sweeps' %blankSweep[0] + ', shuffled'*shuffleBlankSweeps
	#exec("print >> buf, 'Header/sweep checksum: _______ %d / %d' %(headerchksum,sweepchksum)")*DTBOARDINSTALLED
	print >> buf, 'Screen (w,h,d): ______________ (%.1f, %.1f, %.1f) cm, %.2f gamma, %s' %(SCREENWIDTHCM,SCREENHEIGHTCM,SCREENDISTANCECM,GAMMA,EYETEXT[EYE])
	print >> buf, 'Refresh rate: ________________ %s Hz' %REFRESHRATE
	print >> buf, 'Origin: ______________________ (%.1f, %.1f) deg, wrt screen center, center anchor' %(origDeg)
	print >> buf, 'Pre/post-experiment delay: ___ %.3f / %.3f sec' %(preexpSec,postexpSec)
	print >> buf, 'Orientation offset: __________ %.1f deg' %(orioff)
	print >> buf, 'Grid elements (w/h): _________ %d / %d cells' %(ncellswide,ncellshigh)
	print >> buf, 'Grid size (w/h): _____________ %.1f / %.1f deg' %(regionwidthDeg,regionheightDeg)
	print >> buf, 'Bar parameters (constants):'
	exec("print >> buf, '  orientation: _______________ %.1f deg' %ori")*('ori' not in sweepTable)
	exec("print >> buf, '  width: _____________________ %.1f deg' %widthDeg")*('widthDeg' not in sweepTable)
	exec("print >> buf, '  height: ____________________ %.1f deg' %heightDeg")*('heightDeg' not in sweepTable)
	exec("print >> buf, '  brightness (0-1): __________ %.3f' %brightness")*('brightness' not in sweepTable)
        exec("print >> buf, '  colour (R,G,B):   __________ (%.1f, %.1f, %.1f)' %(colour[0],colour[1],colour[2])")*('colour' not in sweepTable)
	exec("print >> buf, '  antialiase (0=no, 1=yes): __ %d' %antialiase")*('antialiase' not in sweepTable)
	exec("print >> buf, '  BG brightness (0-1): _______ %.3f' %bgbrightness")*('bgbrightness' not in sweepTable)
	exec("print >> buf, '  sweep duration: ____________ %.1f msec' %sweeptimeMsec")*('sweeptimeMsec' not in sweepTable)
	exec("print >> buf, '  post-sweep delay (ISI): ____ %.1f msec' %postsweepMsec")*('postsweepMsec' not in sweepTable)

	print >> buf, 'Completed',sweeplisti+1,'of',nsweeps,'sweeps,',nvsyncsdisplayed,'of',nsweepvsyncs,'sweep refreshes' # sweeplisti is 0 based, so add 1
	print >> buf, 'Timing histogram:'
	print >> buf, vsynctimer.log_histogram(),
	if pausehappened or resumehappened:
		print >> buf, 'WARNING: Dimstim was paused at some point'
	if quit:
		print >> buf, 'WARNING: Dimstim was interrupted before completion'
	else:
		print >> buf, 'Dimstim completed successfully'

	# close Eyelink connection - moved to public function in Core.py
        if EYELINKINSTALLED :
            closeEyelink(edfFileName,buf)

	logger = logging.getLogger('VisionEgg') # get the VisionEgg logger
	logger.info('\n' + buf.getvalue()) # log to file, start on a new line after " INFO:"
	print buf.getvalue(), # print to screen
