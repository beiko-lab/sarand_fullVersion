"""
File:		read_files_concatenate.py
Aythor:		Somayeh Kafaie
Date:		March 2021
Purpose:	To read a list of files and concatenate them

To run:
	python read_files_concatenate.py --dir <the directory containing files to be concatenated>
    --out <the output file>
"""
from utils import concatenate_files, extract_files

#files_dir = "/media/Data/PostDoc/Dalhousie/Work/Test2/Experiments/CAMI_H_1/AMR_info/sequences"
#cat_file = "AMR_seqs.fasta"

def read_concatenate_files(files_dir, myfiles):
    myfiles = extract_files(files_dir,'')
    concatenate_files(myfiles, cat_file)

if __name__=="__main__":

	text = 'This code is used to read a list of files and concatenate them'
	parser = argparse.ArgumentParser(description=text)
	parser.add_argument('--dir','-D', type=str, default = '',
		help = 'the address of the directory containing the files to be concatenated')
	parser.add_argument('--out', '-O', type=str, default='',
		help = 'the output file containing all concatenated files')
	args = parser.parse_args()
	if not args.dir:
		logging.error('please enter the path for the directory containing the files')
		import pdb; pdb.set_trace()
		sys.exit()
	if not args.out:
		logging.error('please enter the address of the output file')
		import pdb; pdb.set_trace()
		sys.exit()
    read_concatenate_files(args.dir, args.out)
