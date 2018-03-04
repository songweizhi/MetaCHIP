
list_50 = []
for each in open('/Users/songweizhi/Desktop/50.txt'):
    list_50.append(each.strip())

list_47 = []
for each in open('/Users/songweizhi/Desktop/47.txt'):
    list_47.append(each.strip())

shared = []
uniq_to_50 = []
for each in list_50:
    if each in list_47:
        shared.append(each)
    else:
        uniq_to_50.append(each)

uniq_to_47 = []
for each2 in list_47:
    if each2 not in list_50:
        uniq_to_47.append(each2)


print(sorted(shared))
print(len(shared))
print(sorted(uniq_to_50))
print(len(uniq_to_50))
print(sorted(uniq_to_47))
print(len(uniq_to_47))
