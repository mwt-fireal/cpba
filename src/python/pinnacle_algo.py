import argparse
import re
import os
import sys
import pprint

def parse_block(FILE, is_list=False):

	if is_list == True:
		block = list()
	else:
		block = dict()
	
	sre = r"([^};]+)\s={"
	xre = r"([^};]+)\s=([^;]+)"
	pre = r"([\d\.\-]+)\,([\d\.\-]+)(\,)?"
	ere = r"};"
	
	line = FILE.readline().strip()
	
	while len(line) > 0:

		rr = re.match(ere, line)
		
		if not(rr == None):
			#print "returning with x of %d" % (x+1)
			return block
		
		rr = re.match(sre, line)
		
		if not(rr == None):
			variable = rr.group(1).strip()
			
			is_list = False
			
			if variable.count("[") > 0:
				is_list = True
				variable = variable.replace("[", "").replace("]", "")
			
			if variable.count("#") > 0:
				if not(type(block) == list):
					block = list()
					
			sub_block = parse_block(FILE, is_list)
			
			if type(block) == list:
				block.append(sub_block)
			else:
				#setattr(block, variable, sub_block)
				
				if block.has_key(variable):
					if type(block[variable]) != list:
						ll = list()
						ll.append(block[variable])
						ll.append(sub_block)
						block[variable] = ll
					else:
						block[variable].append(sub_block)
				else:
					block[variable] = sub_block
		else:
			
			rr = re.match(pre, line)
			
			if not(rr == None):
				pl = list()
				
				x = rr.group(1).strip()
				y = rr.group(2).strip()
				
				pl.append(x)
				pl.append(y)
				
				if type(block) == list:
					block.append(pl)
				else:
					print "error, found a point (%s,%s) but block isn't a list" % (x, y)
					sys.exit(1)
				
				line = FILE.readline().strip()
				continue
				
			rr = re.match(xre, line)
			
			if not(rr == None):
				variable = rr.group(1).strip()
				value = rr.group(2).strip()
				
				if variable.startswith("DoseGrid"):
					parts = variable.split(".")
				
					if not(block.has_key("DoseGrid")):
					#if not(hasattr(block, "DoseGrid")):
						dg = dict()
						block['DoseGrid'] = dg
						#setattr(block, 'DoseGrid', dg)
					else:
						dg = block['DoseGrid']
					
					if len(parts) == 3:
						dge_s = parts[1].strip()
						dgv = parts[2].strip()
						
						
						if not(dg.has_key(dge_s)):
						#if not(hasattr(dg, dge_s)):
							dge = dict()
							dg[dge_s] = dge
							#setattr(dg, dge_s, dge)
						else:
							dge = dg[dge_s]
							#dge = getattr(dg, dge_s)

						dge[dgv] = value
						#setattr(dge, dgv, value)
					else:
						dgv = parts[1].strip()
						dg[dgv] = value
						#setattr(dg, dgv, value)
				else:
					if type(block) == list:
						block.append(value)
					else:
						if block.has_key(variable) and type(block[variable]) != list:
							ll = list()
							ll.append(block[variable])
							ll.append(value)
							block[variable] = ll
						else:
							block[variable] = value
			else:
				print "unknown line %s" % line
				sys.exit(1)
		
		line = FILE.readline().strip()
		

def trial_descent(root, trial_num=0):

	found_list = list()

	for key in root.keys():
		path = "%s.%s" % ("TrialList.# \"#%s\"" % trial_num, key)
		
		elem = root[key]
		
		if type(elem) == list or type(elem) == dict:
			descend_path(elem, path, found_list)
		else:
			full_value = "%s = %s" % (path, elem)
			
			if full_value.count(".DoseVolume =") > 0:
				found_list.append((path, elem))
				
			#print full_value
	
	return found_list
	
def descend_path(root, base_path, found_list):

	if type(root) == dict:
		for key in root.keys():
			path = "%s.%s" % (base_path, key)
			elem = root[key]
			
			if type(elem) == list:
				descend_path(elem, base_path, found_list)
			elif type(elem) == dict:
				descend_path(elem, path, found_list)
			else:
				full_value = "%s = %s" % (path, elem)
				#print full_value
				if full_value.count(".DoseVolume =") > 0:
					found_list.append((path, elem))
	elif type(root) == list:
		for key in range(len(root)):
			path = "%s.# \"#%s\"" % (base_path, key)
			elem = root[key]
			
			if type(elem) == list or type(elem) == dict:
				descend_path(elem, path, found_list)
			else:
				full_value = "%s = %s" % (path, elem)
				#print full_value
				if full_value.count(".DoseVolume =") > 0:
					found_list.append((path, elem))
	else:
		print "unknown type from path %s, type %s, root is %s" % (base_path, type(root), root)
		sys.exit(1)
			

def parse_isodose(args):
	isodose_file = "%s/plan.Isodose" % args.root
	
	print isodose_file
	
	sre = r"([^};]+)\s={"
	xre = r"([^};]+)\s=([^;]+)"
	pre = r"([\d\.\-]+)\,([\d\.\-]+)(\,)?"
	ere = r"};"

	if os.path.exists(args.isodose):
		print "overriding \"%s\" with supplied isodose file \"%s\"" % (isodose_file, args.isodose)
		isodose_file = args.isodose

	nn = list()
	
	FILE = open(isodose_file, 'r+')

	iso = dict()
	iso['IsodoseControl'] = dict();
	
	start = iso['IsodoseControl'];
	
	start['NormalizationMode'] = "Absolute"


	line = FILE.readline().strip();
	#skip this line.
	
	line = FILE.readline().strip();
	
	rr = re.match(sre, line)
		
	if not(rr == None):
		variable = rr.group(1).strip()
		block = parse_block(FILE)

		start[variable] = block

		line = FILE.readline().strip()

	return iso	
	
def parse_trial(args):

	trial_file = "%s/plan.Trial" % args.root
	
	print trial_file

	if os.path.exists(args.trial):
		print "overriding \"%s\" with supplied trial file \"%s\"" % (trial_file, args.trial)
		trial_file = args.trial

	FILE = open(trial_file, 'r')	

	#ok, we now can loop with a range variable.
	
	sre = r"([^};]+)\s={"
	ere = r"};"

	line = FILE.readline().strip()
	
	trials = list()
	
	trial = dict()
	
	while len(line) > 0:

		rr = re.match(sre, line)
		
		if not(rr == None):
			variable = rr.group(1).strip()
			block = parse_block(FILE)

			trial[variable] = block
			
			trials.append(trial['Trial'])
			
			line = FILE.readline().strip()
			
			

	#print pprint.pprint(trial)
	#print "The VoxelSize is X = %s, Y = %s, Z = %s" % (trial.DoseGrid.VoxelSize.X, trial.DoseGrid.VoxelSize.Y, trial.DoseGrid.VoxelSize.Z)
	
	#print "we found %s trials in this plan.Trial" % (str(len(trials) + 1))
	
	return trials