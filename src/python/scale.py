#!/usr/bin/python

import pinnacle_algo
import sys
import os
import argparse
import pprint
import re
import os.path


def write_new_isodose_info(scale_factor, isolist, fptr):
	fptr.write("IsodoseControl .LineList .DestroyAllChildren = \"\";\n")
	
	
	for i in range(len(isolist)):
		iso = isolist[i]

		fptr.write("IsodoseControl .LineList .CreateChild = \"\";\n");
		fptr.write("IsodoseControl .LineList .Last .Color = %s;\n" % iso['Color'])
		fptr.write("IsodoseControl .LineList .Last .Display2dOn = %s;\n" % iso['Display2dOn'])
		fptr.write("IsodoseControl .LineList .Last .Display3dMode = %s;\n" % iso['Display3dMode'])
		fptr.write("IsodoseControl .LineList .Last .LineWidth = %s;\n" % iso['LineWidth'])
		fptr.write("IsodoseControl .LineList .Last .IsoValue = %s;\n" % str(round(float(iso['IsoValue']) * scale_factor, 1)))
	
parser = argparse.ArgumentParser(description='Interface with Pinnacle')
parser.add_argument('--root', help='patient data root', default='/pinnacle_patient_expansion/BetaPatients/Institution_323/Mount_0/Patient_566/Plan_0');
parser.add_argument('--trial', help='use this trial file (debugging)', default='plan.Trial');
parser.add_argument('--isodose', help='use this trial file (debugging)', default='plan.Isodose');
parser.add_argument('--trialname', help='trial name', default='Static Field');
parser.add_argument('--scalefactor', help='scale factor', default=2.0);
parser.add_argument('--dryrun', help='scale factor', default=0);

args = parser.parse_args()

trials = pinnacle_algo.parse_trial(args)
iso = pinnacle_algo.parse_isodose(args)

LineList = iso['IsodoseControl']['LineList']['IsodoseLine']

trial_num = 0

for i in range(len(trials)):
	trial = trials[i]
	if trial['Name'].count(args.trialname) > 0:
		trial_num = i
	
print "found trial %s of %s, index %s of %s" % (str(trial_num + 1), str(len(trials) + 1), trial_num, len(trials))
	
trial = trials[trial_num] 
dv_list = pinnacle_algo.trial_descent(trial, trial_num)

pprint.pprint(dv_list)

print 'for every beam, scale the doseVolume, which is found at TrialList.# "#(id)".BeamList.# "#(id)".DoseVolume'

DoseGrid = trial['DoseGrid']

dim_x = DoseGrid['Dimension']['X']
dim_y = DoseGrid['Dimension']['Y']
dim_z = DoseGrid['Dimension']['Z']

#"/home/pinnbeta/cpba/1.1c/src/exe/main --hstep=0 --pinnacle=1 
#--pinnaclepath=$HOME/cpba/temp_files --inputpath=$HOME/cpba/1.1c/InputData --zdepth=23"

outdir = "/tmp/cpba"

pargs = dict()
pargs['x'] = dim_x
pargs['y'] = dim_y
pargs['z'] = dim_z
pargs['factor'] = float(args.scalefactor)
pargs['outputdir'] = outdir

if not(os.path.exists(outdir)):
	os.makedirs(outdir)

reload_script = "%s/SCALE.Script" % outdir

if os.path.exists(reload_script):
	os.remove(reload_script)

fptr = open(reload_script, 'w')



for dv in dv_list:

	path = dv[0]
	bi = dv[1]

	#if not(path.count("BeamList")):
	#	print "skipping %s" % path
	#	continue

	rx = r"([0-9]+)"
	
	rr = re.search(rx, bi)

	if rr == None:
		print "Error: binary index %s did not match extraction regex." % bi
		sys.exit(1)
	
	bin_index = "0%s" % rr.groups(1)
	
	pargs['inputfile'] = "%s/plan.Trial.binary.%s" % (args.root, bin_index)
	pargs['index'] = bin_index
	
	exe = "%s/cpba/1.1c/src/exe/dose_scale" % os.environ['HOME']
	
	for key in pargs.keys():
		exe = "%s --%s=%s" % (exe, key, pargs[key])
	
	print exe
	
	if args.dryrun == 0:
		os.system(exe)
	
	out_file = "%s/dv.%s" % (outdir, bin_index)
	
	pinnacle_line = "%s = \\BOB{B}:%s\\;" % (path, out_file)
	fptr.write("%s\n" % pinnacle_line)

write_new_isodose_info(pargs['factor'], LineList, fptr)

fptr.close()
