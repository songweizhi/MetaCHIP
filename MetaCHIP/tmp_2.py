import os
text_file = '/Users/songweizhi/Desktop/a.txt'
text_file_size = os.stat(text_file).st_size
print(text_file_size)