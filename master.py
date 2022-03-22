total_span = 62
total_width = 30

girder_no = [2, 10, 1]
flange_plate_width = [0.300, 1.200, 0.050]
web_plate_width = [0.500, 3.00, 0.050]
plate_thickness = [12, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50]


def number(var_list):
    var_range = (var_list[1] - var_list[0])
    remainder = var_range % var_list[2]
    var_no = (var_range - remainder) / var_list[2] + 2
    return int(var_no)


class NoOfElements:
    top_flange = number(flange_plate_width) * len(plate_thickness)
    web = number(web_plate_width) * len(plate_thickness)
    bott_flange_1 = number(flange_plate_width) * len(plate_thickness)
    bott_flange_2 = number(flange_plate_width) * len(plate_thickness)
    bott_flange_3 = number(flange_plate_width) * len(plate_thickness)


print(NoOfElements.top_flange * NoOfElements.web * NoOfElements.bott_flange_1 * NoOfElements.bott_flange_2 * NoOfElements.bott_flange_3)






