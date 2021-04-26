
from utils import concatenate_files, extract_files

myfiles = extract_files("/media/Data/PostDoc/Dalhousie/Work/Test2/Experiments/CAMI_H_1/AMR_info/sequences",'')
concatenate_files(myfiles, "AMR_seqs.fasta")
