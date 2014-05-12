#!/usr/bin/python

import pinnacle_algo
import sys
import os
import argparse
import pprint

parser = argparse.ArgumentParser(description='Interface with Pinnacle')
parser.add_argument('--root', help='patient data root', default='/pinnacle_patient_expansion/BetaPatients/Institution_323/Mount_0/Patient_785/Plan_0');
parser.add_argument('--trial', help='use this trial file (debugging)', default='plan.Trial');
parser.add_argument('--trial_number', help='use this trial number', default=0);

args = parser.parse_args()

trials = pinnacle_algo.parse_trial(args)

trial = trials[0]

# trim off the quotes around the value.
electron_eo = eval(trial['BeamList']['Beam']['MachineEnergyName'])

# ControlPoint 
cp = trial['BeamList']['Beam']['CPManager']['CPManagerObject']['ControlPointList'][0]
# Gantry
thetadeg = 180.0 - float(cp['Gantry'])

if thetadeg >= 180.0:
	thetadeg = thetadeg - 180.0

print "ThetaDeg = %s" % thetadeg
print "electron_eo is %s" % electron_eo
map = dict()

map['6 MeV'] = "50"
map['9 MeV'] = "100"
map['12 MeV'] = "150"
map['16 MeV'] = "200"
map['20 MeV'] = "250"


if not(map.has_key(electron_eo)):
	print "unknown mapped energy: %s, defaulting to 100" % electron_eo
	eo = "100"
else:
	eo = map[electron_eo];

print "Using eo of %s" % eo

DoseGrid = trial['DoseGrid']

dx = DoseGrid['VoxelSize']['X']
dz = DoseGrid['VoxelSize']['Y']

print "dx = %s, dz = %s" % (dx, dz)

dim_x = DoseGrid['Dimension']['X']
dim_y = DoseGrid['Dimension']['Y']
dim_z = DoseGrid['Dimension']['Z']

exe = "%s/cpba/1.1c/src/exe/main" % os.environ['HOME']

#"/home/pinnbeta/cpba/1.1c/src/exe/main --hstep=0 --pinnacle=1 
#--pinnaclepath=$HOME/cpba/temp_files --inputpath=$HOME/cpba/1.1c/InputData --zdepth=23"
pargs = dict()

pargs['hstep'] = "0"
pargs['pinnacle'] = "1"
pargs['eo'] = eo
pargs['dx'] = dx
pargs['dz'] = dz
pargs['thetadeg'] = thetadeg
pargs['pinnaclepath'] = "%s/cpba/temp_files" % os.environ['HOME']
pargs['inputpath'] = "%s/cpba/1.1c/InputData" % os.environ['HOME']
pargs['zdepth'] = dim_z

for key in pargs.keys():
	exe = "%s --%s=%s" % (exe, key, pargs[key])

os.system(exe)

