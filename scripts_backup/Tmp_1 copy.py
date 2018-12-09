import time
import multiprocessing


def basic_func(x):
    if x == 0:
        return 'zero'
    elif x % 2 == 0:
        return 'even'
    else:
        return 'odd'


def multiprocessing_func(x):
    time.sleep(2)
    print('{} in a {} number'.format(x, basic_func(x)))


def do_something(each):
    time.sleep(2)
    print('processed: %s' % each)


bin_list = ['bin1.fa', 'bin2.fa', 'bin3.fa', 'bin4.fa', 'bin5.fa', 'bin6.fa', 'bin7.fa', 'bin8.fa', 'bin9.fa', 'bin10.fa']


starttime = time.time()
for each_bin in bin_list:
    time.sleep(2)
    print('processed: %s' % each_bin)
print('That took {} seconds'.format(time.time() - starttime))









starttime = time.time()

processes = []

for each_bin in bin_list:
    p = multiprocessing.Process(target=do_something, args=(each_bin,))
    processes.append(p)


for process in processes:
    process.start()

for process in processes:
    process.join()

print('That took {} seconds'.format(time.time() - starttime))


