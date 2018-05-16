import sys
total = 100000000 + 1
point = total/100 # 100 percent
increment = total/50 # the number of '='
for i in range(total):
    if i % (5 * point) == 0:
        sys.stdout.write("\r[" + "=" * (i/increment) + "_" * ((total - i)/increment) + "]" + str(i/point) + "%")
        sys.stdout.flush()