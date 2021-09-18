import os
import shutil
from time import sleep


def force_create_folder(folder_to_create):

    rm_rd = 0
    while os.path.isdir(folder_to_create) is True:
        shutil.rmtree(folder_to_create, ignore_errors=True)

        if rm_rd >= 10:
            print('Failed in removing %s, program exited!' % folder_to_create)
            exit()

        rm_rd += 1
        sleep(1)

    os.mkdir(folder_to_create)




force_create_folder('/Users/songweizhi/Desktop/hahaha')
