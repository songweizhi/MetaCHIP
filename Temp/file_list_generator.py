import os
import glob

target_files = '/Users/weizhisong/Desktop/*.pdf'
names = [os.path.basename(file_name) for file_name in glob.glob(target_files)]
print(names)
