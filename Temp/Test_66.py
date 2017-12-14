
def get_intersection(list_1, list_2):
    intersect_list = []
    for each in list_1:
        if (each in list_2) and (each not in intersect_list):
            intersect_list.append(each)
    return intersect_list




