from utils.utils import *
import glob

remove_zero_profiles(glob.glob(PROFILES_BLIND_OUTPUT_CSV_DIR + '*'))
