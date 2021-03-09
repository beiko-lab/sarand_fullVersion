
import sys
import os
import copy
import argparse
import enum
import logging
import re


def initialize_logger(output_dir, file_name = 'logfile.log'):
	"""
	"""
	#logging_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'logs')
	logging_dir = os.path.join(output_dir,'logs')
	try:
		os.makedirs(logging_dir)
	except OSError:
		pass
	log_file_path = os.path.join(logging_dir, file_name)
	logging.basicConfig(level=logging.DEBUG)
	#logging.basicConfig()
	log_formatter = logging.Formatter(
		"%(asctime)s [%(levelname)-5.5s]	%(message)s"
	)
	root_logger = logging.getLogger()
	file_handler = logging.FileHandler(log_file_path)
	file_handler.setFormatter(log_formatter)
	root_logger.addHandler(file_handler)
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setFormatter(log_formatter)
	root_logger.addHandler(console_handler)

def check_reads(v):
	if isinstance(v, str):
		return v
	elif isinstance(v, list):
		if len(v)!=2:
			return False
		for item in v:
			if not isinstance(v, str):
				return False
		return v
	else:
		return False

def str2bool(v):
	"""
	To convert a string to a boolean value
	Parameter:
		v: the input string
	Return:
		the converted boolean value
	"""
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')

def print_params(params):
	"""
	"""
	mylog = '\n'.join(param for param in params)
	logging.info(mylog)

def verify_file_existence(myfile, param_file, message):
	"""
	This code is used to check if the required file exits either through running
	previous steps of the pipeline or passed by user (or in the param list)
	Parameters:
		myfile:	the file was supposed to be generated via previous steps of the pipeline
		param_file:	the file set by the user
		message:	error message in case no instance of the required file exists
	Return:
		the valid instance of the parameter to work with
	"""
	if myfile!="":
		return myfile
	elif param_file !="":
		if isinstance(param_file, str) and os.path.isfile(param_file):
			return param_file
		elif isinstance(param_file, list):
			error = False
			for file in param_file:
				if not isinstance(file, str) or not os.path.isfile(file):
					error = True
			if not error:
				return True
	logging.error("ERROR: "+message)
	#print("ERROR: "+message)
	sys.exit()

def retrieve_AMR(file_path):
	"""
	To read the AMR gene from the text file.
	Parameters:
		file_path:	the address of the file containing the AMR gene
	Return:
		the sequence of the AMR gene in lower case
	"""
	with open(file_path) as fp:
		for i, line in enumerate(fp):
			#skip comment line
			if line.startswith('>'):
				continue
			return line

def reverse_sign(sign):
	"""
	To reverse the sign (+/-)
	Parameetr:
		sign: either + or -
	Return:
		the reverse sign
	"""
	if sign=='-':
		return '+'
	elif sign=='+':
		return '-'
	else:
		logging.error("ERROR: ivalid sign!")
		#print("ERROR: ivalid sign!")
		sys.exit()

def find_node_orient(node):
	"""
	To remove specific characters and return the last character of what remains
	as the orient of the node
	"""
	return re.sub('[]}]', '', node)[-1]

def find_node_name(node):
	"""
	To remove specific characters and return the rest except the last character
	as the node name
	"""
	return re.sub('[]{}[]','', node)[:-1]

def find_node_name_orient(node):
	"""
	To remove specific characters and return the rest as the name+orient of the node
	"""
	return re.sub('[]{}[]','', node)

def exist_in_path(path, mynode):
	"""
	To check if a given node exists in the path
	Parameters:
		path:	a list of nodes
		mynde:	the node to check if it is in the path
	Return:
		True if mynode is in the path; False otherwise
	"""
	for i, node in enumerate(path):
		if find_node_name_orient(node)==mynode:
			return i
	return -1

def extract_files(gfiles, message):
	"""
	To extract file(s) address from an object
	# if gfiles is a list and the first item is a file address (it would be more
	# accurate to check this for all items) return gfiles
	# else if gfiles is a file address return [gfiles] as a list with one item
	# else if gfiles is a directory address return the list of files in it
	Parameters:
		gfiles:		a string or list of strings
		message:	an error message in case that no file was extracted successfully
	Return:
		the list of file(s) address
	"""
	if isinstance(gfiles, list):
		#check if these are files (not directories)
		if os.path.isfile(gfiles[0]):
			return gfiles
		else:
			logging.error(message)
			#print(message)
			sys.exit()
		#elif os.path.isdir(gfiles[0])
	elif os.path.isfile(gfiles):
		return [gfiles]
	elif os.path.isdir(gfiles):
		myfiles = [os.path.join(gfiles, f) for f in os.listdir(gfiles) \
							if os.path.isfile(os.path.join(gfiles, f))]
		return myfiles
	elif message!='':
		logging.error(message)
		#print(message)
		sys.exit()
