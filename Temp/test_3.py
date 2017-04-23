
input_list = [['None'], ['None'], ['None']]


def uniq_list(input_list):
    output_list = []
    for each_element in input_list:
        if each_element not in output_list:
            output_list.append(each_element)
    return output_list


output_list = uniq_list(input_list)
print(output_list)