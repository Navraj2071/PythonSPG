

def react(no_of_gir, spacing, cantilever, load, loadoffset):
    fem = [[0, 0]]
    dem = [[0, 0]]
    tem = [[0, 0]]
    factor = [[0, 0]]
    reactions = [0]
    d = [0]
    for i in range(no_of_gir):
        fem.insert(i, [0, 0])
        dem.insert(i, [0, 0])
        tem.insert(i, [0, 0])
        factor.insert(i, [0, 0])
        reactions.insert(i, 0)
        d.insert(i, 0)
    if no_of_gir <= 3:
        p = 0.5
        q = 0.5
    else:
        p = 3 / 7
        q = 4 / 7

    if loadoffset > cantilever * 2 + (no_of_gir - 1) * spacing:
        return 0
    else:
        pass

    loadlocation = 0
    if loadoffset < cantilever:
        loadlocation = 0
    else:
        i = 0
        while loadoffset > cantilever + spacing * i:
            loadlocation += 1
            i += 1
    if loadlocation == 0:
        fem[0][1] = load * (cantilever - loadoffset)
    elif loadlocation == no_of_gir:
        fem[no_of_gir][0] = -1 * load * (loadoffset - cantilever - (no_of_gir - 1) * spacing)
    else:
        for i in range(1, no_of_gir):
            if loadlocation == i:
                a = loadoffset - cantilever - (i - 1) * spacing
                b = cantilever + i * spacing - loadoffset
                fem[i][0] = -1 * load * a * b ** 2 / spacing ** 2
                fem[i][1] = load * a ** 2 * b / spacing ** 2
            else:
                pass

    dem[0][0] = fem[0][0]
    dem[0][1] = fem[0][1]
    dem[no_of_gir][0] = fem[no_of_gir][0]
    dem[no_of_gir][1] = fem[no_of_gir][1]
    dem[1][0] = dem[0][1] * -1
    dem[no_of_gir - 1][1] = dem[no_of_gir][0] * -1
    fem[1][1] -= (fem[1][0] / 2 + fem[0][1] / 2)
    fem[no_of_gir - 1][0] -= (fem[no_of_gir - 1][1] / 2 + fem[no_of_gir][0] / 2)

    for i in range(1, no_of_gir):
        factor[i] = [0.5, 0.5]
    factor[1][1] = p
    factor[2][0] = q
    factor[no_of_gir - 2][1] = q
    factor[no_of_gir - 1][0] = p
    mom_diff_limit = 0.001
    is_iterate = True
    iteration_no = 1
    while is_iterate:
        mom_diff = 0
        for i in range(1, no_of_gir):
            fem[0][0] = 0
            fem[0][1] = 0
            fem[no_of_gir][0] = 0
            fem[no_of_gir][1] = 0
            fem[1][0] = 0
            fem[no_of_gir - 1][1] = 0

            tem[i][0] = fem[i][0] - (fem[i - 1][1] + fem[i][0]) * factor[i][0] - \
                        (fem[i][1] + fem[i + 1][0]) * factor[i][1] * 0.5
            tem[i][1] = fem[i][1] - (fem[i][1] + fem[i + 1][0]) * factor[i][1] - \
                        (fem[i - 1][1] + fem[i][0]) * factor[i][0] * 0.5
        for i in range(1, no_of_gir):
            fem[i][0] = tem[i][0]
            fem[i][1] = tem[i][1]
        for i in range(1, no_of_gir):
            if abs(abs(fem[i][1]) - abs(fem[i + 1][0])) > mom_diff:
                mom_diff = abs(abs(fem[i][1]) - abs(fem[i + 1][0]))
        if mom_diff > mom_diff_limit:
            is_iterate = True
        else:
            is_iterate = False
        if iteration_no > 100:
            is_iterate = False
        iteration_no += 1

    dem[1][1] = fem[1][1]
    dem[no_of_gir - 1][0] = fem[no_of_gir - 1][0]

    for i in range(2, no_of_gir - 1):
        dem[i][0] = fem[i][0]
        dem[i][1] = fem[i][1]

    for i in range(1, no_of_gir):
        girderoffset = cantilever + i * spacing
        loaddistance = girderoffset - loadoffset
        if loadoffset > girderoffset:
            loaddistance = 0
        else:
            pass
        reaction_sum = 0
        for n in range(i - 1):
            d[n] = spacing * (i - n)
            reaction_sum += d[n] * reactions[n]
        reactions[i - 1] = (-1 * dem[i][1] + load * loaddistance - reaction_sum) / spacing
    reaction_sum = sum(reactions)
    reactions[no_of_gir - 1] = load - reaction_sum
    return reactions


# srty = react(10, 3, 1.5, 50.0, 18.69)
# print(srty)
