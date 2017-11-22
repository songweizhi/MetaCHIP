

def turn_to_percentage(number_list):
    number_list_percent = []
    for each_element in number_list:
        each_element_percent = float("{0:.2f}".format(each_element / sum(number_list)))
        number_list_percent.append(each_element_percent)
    return number_list_percent


number_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

print(turn_to_percentage(number_list))
