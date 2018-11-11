# PATHS
INPUT_DSSP_FOLDER = 'data/input/dssp/'
INPUT_FASTA_FOLDER = 'data/input/fasta/'
INPUT_CV_DIR = 'data/input/svm/cv/'
TRAIN_FASTA_FILE_NAME = "data/blind_test/train.fasta"
ROW_BLIND_FASTA_FILE_NAME = 'data/input/blind_test/row_blind.fasta'
BLIND_FASTA_FILE_NAME = 'data/blind_test/blind.fasta'
BLIND_DSSP_FILE_NAME = 'data/blind_test/blind.dssp'
BLIND_DSSP_DIR = 'data/blind_test/out/dssp/'
BLIND_FASTA_DIR = 'data/blind_test/out/fasta/'
BLASTCLUST_OUT_FILE_NAME = 'data/blind_test/blastclust.out'
BLASTP_OUT_FILE_NAME = 'data/blind_test/blastp_out.tab'
PROFILES_DB = "data/input/profiles/uniprot_sprot.fastadb"
PROFILES_OUTPUT_PSSM_DIR = "data/profiles/out/pssm/"
PROFILES_OUTPUT_BLAST_DIR = "data/profiles/out/blast/"
PROFILES_OUTPUT_CSV_DIR = "data/profiles/out/csv_profiles/"
PROFILES_BLIND_OUTPUT_PSSM_DIR = "data/profiles_blind/out/pssm/"
PROFILES_BLIND_OUTPUT_BLAST_DIR = "data/profiles_blind/out/blast/"
PROFILES_BLIND_OUTPUT_CSV_DIR = "data/profiles_blind/out/csv_profiles/"
PROFILES_BLIND_INPUT_FASTA_FILE = "data/blind_test/blind.fasta"
GOR_OUT_DIR = "data/gor/"
GOR_TRAINED_FILE = GOR_OUT_DIR + "gor_trained.csv"
SVM_DIR = "data/svm/"
SVM_TESTS_DIR = "data/svm/tests/"
SVM_BLIND_VECTORS_DIR = "data/svm/blind_vectors/"

SCRIPT_MAKEBLASTDB = 'scripts/makeblastdb.sh'
SCRIPT_BLASTCLUST = 'scripts/blastclust.sh'
SCRIPT_BLASTP = 'scripts/blastp.sh'
SCRIPT_PSIBLAST = 'scripts/psiblast.sh'

# URLS
PDB_URL = 'http://files.rcsb.org/download/'

# STRING CONSTANTS
TMP = '_tmp'

# common order od secondary structure keys
SS_ORDER = ['h', 'e', 'c', 't']
# common order of amino acids
AA_ORDER = ['a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v']
AA_ORDER_STR = 'arndcqeghilkmfpstwyv'
# second structure mapping
SS_MAP = {'H': 'H',
          'G': 'H',
          'I': 'H',
          'B': 'E',
          'E': 'E',
          'T': 'C',
          'S': 'C',
          '-': 'C'}
# window size for GOR
window = 17
wind_l = 0 - window // 2
wind_r = window // 2

C_GAMMA_SET = [(1, 1), (2, 0.8), (2.6, 0.5), (1.8, 0.2)]
