import codedata
import master
import trans_dist
import math


def design_run(top_flange_width):
    # ##############################################INPUTS################
    class Deck:
        width = 30
        length = 62
        thickness = 0.22
        haunch_thickness = 0.15
        haunch_width = 0.15
        long_cantilever = 0.5
        fck = 40 * 100
        density = 2.5
        green_density = 2.6
        cb_left = 0.5
        cb_right = 0.5


    class Girder:
        no = 10
        spacing = 3
        overhang = 0.5
        fy = 330 * 100
        fu = 490 * 100
        density = 7.85
        elastic_modulus = 200000 * 100
        rigidity_modulus = 77000 * 100
        mu = 0.3


    class Parameters:
        m_long_term = 15
        m_short_term = 75


    components = ['top_flange', 'web', 'bott_flange_1', 'bott_flange_2', 'bott_flange_3']
    height = {
        components[0]: 0.02,
        components[1]: 1.5,
        components[2]: 0.025,
        components[3]: 0.05,
        components[4]: 0.032
    }

    width = {
        components[0]: 1.45,
        components[1]: 0.044,
        components[2]: 0.95,
        components[3]: 0.7,
        components[4]: 1.0
    }


    class Loads:
        shuttering = 1.5
        const_ll = 0.3
        surfacing = 0.2
        sidl = 1
        sidl_location = [0, Deck.width]
        basic_wind_speed = 56
        super_height = 100
        cb_height = 1
        terrain = 'Plain'
        seis_coeff = 0.72
        sv_load = True


    cantilever = (Deck.width - (Girder.no - 1) * Girder.spacing) / 2


    # ############################## PROPERTIES
    def elastic_properties(height, width, components):
        cg = [0]
        area = [0]
        second_moment_of_inertia_zz = [0]
        second_moment_of_inertia_yy = [0]
        sf_flange = [0]
        ecc_from_topflange = [0]

        area_moment = 0
        area_of_section = 0
        dna = 0
        second_moment_of_inertia_zz_of_section = 0
        second_moment_of_inertia_yy_of_section = 0
        shear_center = 0

        for i in range(len(components)):
            # Populate list of cg of every member
            if i == 0:
                cg[i] = height[components[i]] / 2
            else:
                cg.insert(i, cg[i - 1] + height[components[i]] / 2 + height[components[i - 1]] / 2)
            # Populate list of shear force in flange of every member
            if components[i] != 'web':
                sf_flange.insert(i, height[components[i]] * width[components[i]] ** 3)
            else:
                sf_flange.insert(i, 0)

            ecc_from_topflange.insert(i, cg[i] - cg[0])
            # Populate list of area of cross section of every member
            area.insert(i, height[components[i]] * width[components[i]])
            # Populate list of area of local moment of inertia of every member
            second_moment_of_inertia_zz.insert(i, width[components[i]] * height[components[i]] ** 3 / 12)
            second_moment_of_inertia_yy.insert(i, width[components[i]] ** 3 / 12 * height[components[i]])

            # Total moment of area of every section
            area_moment += area[i] * cg[i]
            # Total area of Girder Section
            area_of_section += area[i]
            # Eccentricity of shear center x 12 Iy
            shear_center += sf_flange[i] * ecc_from_topflange[i]

        dna = area_moment / area_of_section
        girder_height = 0
        for i in range(len(components)):
            girder_height += height[components[i]]
            second_moment_of_inertia_zz[i] += area[i] * (cg[i] - dna) ** 2
            second_moment_of_inertia_zz_of_section += second_moment_of_inertia_zz[i]
            second_moment_of_inertia_yy_of_section += second_moment_of_inertia_yy[i]

        shear_center = shear_center / 12 / second_moment_of_inertia_yy_of_section
        output = {
            'area of cross section': area_of_section,
            'moment of inertia zz': second_moment_of_inertia_zz_of_section,
            'moment of inertia yy': second_moment_of_inertia_yy_of_section,
            'shear center': shear_center,
            'depth of na': dna,
            'section modulus': second_moment_of_inertia_zz_of_section / dna,
            'yg': girder_height / 2 - (shear_center + height[components[0]] / 2)
        }
        return output


    def plastic_properties(height, width, components, comp_uls, tens_uls):
        dpoint = [0]
        comp_area = [0]
        tens_area = [0]
        comp_leverarm = [0]
        tens_leverarm = [0]
        section_modulus = 0
        girder_height = 0

        for i in range(len(components) - 1):
            comp_area.insert(i, 0)
            tens_area.insert(i, 0)
            comp_leverarm.insert(i, 0)
            tens_leverarm.insert(i, 0)

        for i in range(len(components)):
            girder_height += height[components[i]]
            dpoint.insert(i, girder_height)
        net_area = 10
        dpna = girder_height / 2
        max_itr = 100
        itr_count = 0
        incr = dpna / 2
        while 0.00001 < net_area or net_area < -0.00001:
            section_modulus = 0
            for i in range(len(components)):
                if i == 0:
                    if height[components[i]] > dpna:
                        area_in_compression = dpna * width[components[i]]
                        area_in_tension = (height[components[i]] - dpna) * width[components[i]]
                        comp_leverarm[i] = dpna / 2
                        tens_leverarm[i] = (height[components[i]] - dpna) / 2
                    else:
                        area_in_compression = height[components[i]] * width[components[i]]
                        area_in_tension = 0
                        comp_leverarm[i] = dpna - height[components[i]] / 2
                        tens_leverarm[i] = 0
                else:
                    if dpoint[i] > dpna:
                        area_in_compression = max(0, (dpna - dpoint[i - 1])) * width[components[i]]
                        area_in_tension = (dpoint[i] - max(dpoint[i - 1], dpna)) * width[components[i]]
                        comp_leverarm[i] = (dpna - dpoint[i - 1]) / 2
                        tens_leverarm[i] = (dpoint[i] - max(dpoint[i - 1], dpna)) / 2 + max(dpoint[i - 1] - dpna, 0)
                    else:
                        area_in_compression = height[components[i]] * width[components[i]]
                        area_in_tension = 0
                        comp_leverarm[i] = dpna - dpoint[i] + height[components[i]] / 2
                        tens_leverarm[i] = 0
                comp_area[i] = area_in_compression * comp_uls[i]
                tens_area[i] = area_in_tension * tens_uls[i]

                section_modulus += comp_area[i] * comp_leverarm[i] + tens_area[i] * tens_leverarm[i]
            net_area = sum(comp_area) - sum(tens_area)
            if net_area > 0:
                dpna = dpna - incr
            else:
                dpna = dpna + incr
            incr = incr / 2

            if itr_count > max_itr:
                net_area = 0
            itr_count = itr_count + 1
        output = {
            'depth of plastic na': dpna,
            'section modulus': section_modulus,
        }
        return output


    def section_classification(fy, width, thickness):
        epsilon = (250 / fy) ** 0.5
        limiting_ratio = [8.4 * epsilon, 9.4 * epsilon, 13.6 * epsilon, 100]
        classification = ['Plastic', 'Compact', 'Semi-Compact', 'Slender']
        i = 0
        while width / thickness > limiting_ratio[i]:
            i += 1
        return classification[i]


    def girder_bending_strength(height, width, components, unsupp_length, fy, elastic_modulus_steel,
                                rigidity_modulus_steel, elastic_output, plastic_output):
        elastic_section_modulus = elastic_output['section modulus']
        plastic_section_modulus = plastic_output['section modulus']
        section_class = section_classification(fy, (width[components[0]] - width[components[1]]) / 2, height[components[0]])
        if section_class == 'Semi compact':
            beta_b = elastic_section_modulus / plastic_section_modulus
        else:
            beta_b = 1.0
        factor = {
            'K': 1,
            'Kw': 1,
            'c1': 1.132,
            'c2': 0.459,
            'c3': 0.525
        }
        yg = elastic_output['yg']
        moi_comp_flange = width[components[0]] * height[components[0]] ** 3 / 12
        moi_tens_flange = width[components[-1]] * height[components[-1]] ** 3 / 12 + \
                          width[components[-2]] * height[components[-2]] ** 3 / 12 + \
                          width[components[-3]] * height[components[-3]] ** 3 / 12
        beta_f = moi_comp_flange / (moi_comp_flange + moi_tens_flange)
        hy = height[components[0]] / 2 + height[components[1]] + (height[components[2]] + height[components[3]] +
                                                                  height[components[4]]) / 2
        if beta_f > 0.5:
            yj = 0.8 * (2 * beta_f - 1) * hy / 2
        else:
            yj = (2 * beta_f - 1) * hy / 2
        warping_constant = (1 - beta_f) * beta_f * elastic_output['moment of inertia yy'] * hy ** 2
        torsion_constant = 0
        for i in components:
            torsion_constant += max(width[i], height[i]) * min(width[i], height[i]) ** 3 / 3
        effective_length = 1.2 * unsupp_length
        mcr1 = factor['c1'] * math.pi ** 2 * elastic_modulus_steel * elastic_output['moment of inertia yy'] / (
                effective_length ** 2)
        mcr2 = (factor['K'] / factor['Kw']) ** 2 * warping_constant / elastic_output['moment of inertia yy']
        mcr3 = rigidity_modulus_steel * torsion_constant * effective_length ** 2 / (
                elastic_modulus_steel * elastic_output['moment of inertia yy'] * math.pi ** 2)
        mcr4 = (factor['c2'] * yg - factor['c3'] * yj)
        mcr = mcr1 * ((mcr2 + mcr3 + mcr4 ** 2) ** 0.5 - mcr4)
        lambda_lt = min((beta_b * plastic_section_modulus * fy / mcr) ** 0.5,
                        (1.2 * elastic_section_modulus * fy / mcr) ** 0.5)
        alpha_lt = 0.49
        phi_lt = 0.5 * (1 + alpha_lt * (lambda_lt - 0.2) + lambda_lt ** 2)
        chi_lt = min(1 / (phi_lt + (phi_lt ** 2 - lambda_lt ** 2) ** 0.5), 1.0)
        fbd = chi_lt * fy / 1.1
        bending_strength = beta_b * plastic_section_modulus * fbd
        return bending_strength


    def girder_shear_strength(height, width, stiff_spacing, fy, elastic_modulus_steel, mu):
        kv = 5.35 + 4 / (stiff_spacing / height['web']) ** 2
        tao_cr = kv * math.pi ** 2 * elastic_modulus_steel / (12 * (1 - mu ** 2) * (height['web'] / width['web']) ** 2)
        lambda_w = (fy / (tao_cr * 3 ** 0.5)) ** 0.5
        if lambda_w <= 0.8:
            tao_b = fy / 3 ** 0.5
        elif lambda_w < 1.2:
            tao_b = (1 - 0.8 * (lambda_w - 0.8) * fy / 3 ** 0.5)
        else:
            tao_b = fy / (lambda_w * 3 ** 0.5)
        shear_area = width['web'] * height['web']
        shear_strength_post_critical = shear_area * tao_b
        shear_strength_plastic = shear_area * fy / 3 ** 0.5
        return min(shear_strength_plastic, shear_strength_post_critical)


    stress_uls = [1.0] * len(height)
    girder_depth = 0
    for i in range(len(components)):
        stress_uls.insert(i, 1)
        girder_depth += height[components[i]]

    elastic_properties_output = elastic_properties(height, width, components)
    plastic_properties_output = plastic_properties(height, width, components, stress_uls,
                                                   stress_uls)
    eff_span = Deck.length - Deck.long_cantilever * 2 - Girder.overhang * 2
    stiff_spacing = 4.0


    class PlateGirder:
        area = elastic_properties_output['area of cross section']
        depth = girder_depth
        bending_strength_stage1 = girder_bending_strength(height, width, components,
                                                          eff_span, Girder.fy, Girder.elastic_modulus,
                                                          Girder.rigidity_modulus, elastic_properties_output,
                                                          plastic_properties_output)
        bending_strength_stage2 = girder_bending_strength(height, width, components,
                                                          stiff_spacing, Girder.fy, Girder.elastic_modulus,
                                                          Girder.rigidity_modulus, elastic_properties_output,
                                                          plastic_properties_output)
        shear_strength = girder_shear_strength(height, width, stiff_spacing, Girder.fy,
                                               Girder.elastic_modulus, Girder.mu)


    # ###################### COMPOSITE PROPERTIES
    eff_width_inn_gir = min(eff_span / 4, Girder.spacing)
    eff_width_out_gir = min(Girder.spacing / 2, eff_span / 8) + min(cantilever, eff_span / 8)

    comp_components = ['deck', 'haunch', components[0], components[1], components[2], components[3], components[4]]
    comp_height = {
        comp_components[0]: Deck.thickness,
        comp_components[1]: Deck.haunch_thickness,
        comp_components[2]: height[comp_components[2]],
        comp_components[3]: height[comp_components[3]],
        comp_components[4]: height[comp_components[4]],
        comp_components[5]: height[comp_components[5]],
        comp_components[6]: height[comp_components[6]]
    }

    outer_width = {
        comp_components[0]: eff_width_out_gir,
        comp_components[1]: width[comp_components[2]],
        comp_components[2]: width[comp_components[2]],
        comp_components[3]: width[comp_components[3]],
        comp_components[4]: width[comp_components[4]],
        comp_components[5]: width[comp_components[5]],
        comp_components[6]: width[comp_components[6]]
    }

    inner_width = {
        comp_components[0]: eff_width_inn_gir,
        comp_components[1]: width[comp_components[2]],
        comp_components[2]: width[comp_components[2]],
        comp_components[3]: width[comp_components[3]],
        comp_components[4]: width[comp_components[4]],
        comp_components[5]: width[comp_components[5]],
        comp_components[6]: width[comp_components[6]]
    }

    comp_uls = [Deck.fck * 0.446, Deck.fck * 0.446, Girder.fy / 1.1, Girder.fy / 1.1, Girder.fy / 1.1, Girder.fy / 1.1,
                Girder.fy / 1.1]

    tens_uls = [0, 0, Girder.fy / 1.1, Girder.fy / 1.1, Girder.fy / 1.1,
                Girder.fy / 1.1, Girder.fy / 1.1]

    comp_girder_depth = 0
    for i in comp_components:
        comp_girder_depth += comp_height[i]

    outer_plastic_properties = plastic_properties(comp_height, outer_width, comp_components, comp_uls, tens_uls)
    inner_plastic_properties = plastic_properties(comp_height, inner_width, comp_components, comp_uls, tens_uls)

    outer_dpna = outer_plastic_properties['depth of plastic na']
    outer_mor = outer_plastic_properties['section modulus']

    inner_dpna = inner_plastic_properties['depth of plastic na']
    inner_mor = inner_plastic_properties['section modulus']


    class CompGirder:
        depth = comp_girder_depth
        outer_bending_strength = outer_mor
        inner_bending_strength = inner_mor


    # ################################# LOADS
    sidl_list = [0] * len(Loads.sidl_location)
    for i in range(len(Loads.sidl_location)):
        sidl_list[i] = trans_dist.react(Girder.no, Girder.spacing, cantilever, Loads.sidl, Loads.sidl_location[i])

    sidl_inner_per_girder = [0] * (Girder.no - 2)
    for i in range(Girder.no - 2):
        for j in range(len(Loads.sidl_location)):
            sidl_inner_per_girder[i] += sidl_list[j][i + 1]

    sidl_outer_per_girder = [0] * 2
    for j in range(len(Loads.sidl_location)):
        sidl_outer_per_girder[0] += sidl_list[j][0]
        sidl_outer_per_girder[1] += sidl_list[j][Girder.no - 1]


    def wind_vertical(basic_wind_speed, super_height, terrain):
        i = 0
        while super_height > codedata.super_height[i]:
            i += 1
        if terrain == 'Plain':
            mean_wind_speed = ((codedata.mean_wind_speed_plain[i + 1] - codedata.mean_wind_speed_plain[i]) / (
                    codedata.super_height[i + 1] - codedata.super_height[i]) * (super_height - codedata.super_height[i]) +
                               codedata.mean_wind_speed_plain[i]) / 33 * basic_wind_speed
        else:
            mean_wind_speed = ((codedata.mean_wind_speed_obs[i + 1] - codedata.mean_wind_speed_obs[i]) / (
                    codedata.super_height[i + 1] - codedata.super_height[i]) * (super_height - codedata.super_height[i]) +
                               codedata.mean_wind_speed_obs[i]) / 33 * basic_wind_speed
        gust_factor = 2
        coeff_lift = 0.75
        hor_wind_press = 0.6 * mean_wind_speed ** 2 / 10000
        load_vertical = hor_wind_press * gust_factor * coeff_lift
        return load_vertical


    class Udl:
        girder_weight = PlateGirder.area * Girder.density * 1.35
        deck_weight_outer = (eff_width_out_gir * Deck.thickness +
                             Deck.haunch_thickness * width[components[0]] +
                             Deck.haunch_thickness * Deck.haunch_width) * Deck.density
        deck_weight_inner = (eff_width_inn_gir * Deck.thickness +
                             Deck.haunch_thickness * width[components[0]] +
                             Deck.haunch_thickness * Deck.haunch_width) * Deck.density
        deck_weight_outer_green = deck_weight_outer / Deck.density * Deck.green_density
        deck_weight_inner_green = deck_weight_inner / Deck.density * Deck.green_density
        shutt_outer = eff_width_out_gir * Loads.shuttering
        shutt_inner = eff_width_inn_gir * Loads.shuttering
        const_ll_outer = eff_width_out_gir * Loads.const_ll
        const_ll_inner = eff_width_inn_gir * Loads.const_ll
        surf_outer = eff_width_out_gir * Loads.surfacing
        surf_inner = eff_width_inn_gir * Loads.surfacing
        sidl_inner = max(sidl_inner_per_girder)
        sidl_outer = max(sidl_outer_per_girder)
        wind_const_outer = wind_vertical(Loads.basic_wind_speed, Loads.super_height, Loads.terrain) * eff_width_out_gir * 0.7
        wind_const_inner = wind_vertical(Loads.basic_wind_speed, Loads.super_height, Loads.terrain) * eff_width_inn_gir * 0.7
        wind_noll_outer = wind_const_outer / 0.7
        wind_noll_inner = wind_const_inner / 0.7
        wind_ll_outer = wind_vertical(min(Loads.basic_wind_speed, 36), Loads.super_height,
                                      Loads.terrain) * eff_width_out_gir
        wind_ll_inner = wind_vertical(min(Loads.basic_wind_speed, 36), Loads.super_height,
                                      Loads.terrain) * eff_width_inn_gir
        seis_noll_outer = (girder_weight + deck_weight_outer + surf_outer + sidl_outer) * Loads.seis_coeff
        seis_noll_inner = (girder_weight + deck_weight_inner + surf_inner + sidl_inner) * Loads.seis_coeff


    cb_left = Deck.cb_left
    cb_right = Deck.cb_right
    clear_carriageway = Deck.width - cb_left - cb_right
    i = 0
    while clear_carriageway >= codedata.carriageway[i]:
        i += 1
    no_lane = codedata.lanes[i - 1]
    max_clssA = no_lane
    max_70r = (no_lane - (no_lane % 2)) / 2


    def ll_react(lanes, load_width, min_pos, max_pos, lane_dis, no_of_girder, spacing, cantilever):
        increment = (max_pos - min_pos) / 10
        load_pos = [0]
        for i in range(11):
            t = [min_pos + increment * i, min_pos + load_width + increment * i]
            for j in range(int(lanes * 2 - 2)):
                t.insert(j + 2, (t[j] + lane_dis))
            load_pos.insert(i, t)
        reaction = [[0] * no_of_girder for _ in range(11)]
        load = 50.0

        for i in range(11):
            for k in range(int(lanes * 2)):
                reactions = trans_dist.react(no_of_girder, spacing, cantilever, load, load_pos[i][k])
                for j in range(no_of_girder):
                    reaction[i][j] += reactions[j]
        output = []
        for i in range(no_of_girder):
            maxima = 0
            for k in range(11):
                if reaction[k][i] > maxima:
                    maxima = reaction[k][i]
                else:
                    pass
            output.append(maxima)
        return output


    lcase = [0] * int(2 + max_70r)
    for i in range(int(1 + max_70r)):
        lcase[i] = [max_clssA - 2 * i, i]
    lcase[int(1 + max_70r)] = [0, max_70r]


    def sv_load():
        sv_width = 1.8
        case1_pos_t1 = cb_left + clear_carriageway / 2 - sv_width / 2
        case1_pos_t2 = case1_pos_t1 + sv_width
        case2_pos_t1 = case1_pos_t1 + 0.3
        case2_pos_t2 = case2_pos_t1 + sv_width
        react_case1_t1 = trans_dist.react(Girder.no, Girder.spacing, cantilever, 50, case1_pos_t1)
        react_case1_t2 = trans_dist.react(Girder.no, Girder.spacing, cantilever, 50, case1_pos_t2)
        react_case2_t1 = trans_dist.react(Girder.no, Girder.spacing, cantilever, 50, case2_pos_t1)
        react_case2_t2 = trans_dist.react(Girder.no, Girder.spacing, cantilever, 50, case2_pos_t2)

        react_max = [0] * Girder.no
        for i in range(Girder.no):
            react_max[i] = max(react_case2_t1[i], react_case2_t2[i], react_case1_t1[i], react_case1_t2[i])
        return react_max


    def case1_clA(lcase_list):
        react_max = [0] * Girder.no
        for i in range(len(lcase_list)):
            clssA_lanes = lcase_list[i][0]
            clss70r_lanes = lcase_list[i][1]
            clssA_width = 1.8
            clssA_react = [0] * Girder.no
            if clssA_lanes == 0:
                pass
            else:
                clssA_min_pos = cb_left + 0.15 + 0.25
                if clss70r_lanes == 0:
                    clssA_max_pos = Deck.width - cb_right - 0.15 - 0.25 - 1.8 - (clssA_lanes - 1) * 3.5
                else:
                    clssA_max_pos = Deck.width - cb_right - 7.25 - (clss70r_lanes - 1) * 7 - 1.2 - 0.25 - 1.8 - \
                                    (clssA_lanes - 1) * 3.5
                clssA_react = ll_react(clssA_lanes, clssA_width, clssA_min_pos, clssA_max_pos, 3.5, Girder.no,
                                       Girder.spacing, cantilever)
                for k in range(Girder.no):
                    if clssA_react[k] > react_max[k]:
                        react_max[k] = clssA_react[k]
        return react_max


    def case1_70R(lcase_list):
        react_max = [0] * Girder.no
        for i in range(len(lcase_list)):
            clss70r_lanes = lcase_list[i][1]
            clss70r_width = 1.93
            clss70r_react = [0] * Girder.no
            if clss70r_lanes == 0:
                pass
            else:
                clss70r_react_per_lane = [0] * Girder.no
                for j in range(1, int(clss70r_lanes) + 1):
                    clss70r_min_pos = Deck.width - cb_right - 7.25 - (j - 1) * 7 + 0.43
                    if j == 1:
                        clss70r_max_pos = Deck.width - cb_right - 1.2 - 0.43 - clss70r_width
                    elif j > 1:
                        clss70r_max_pos = Deck.width - cb_right - 7.25 - (j - 2) * 7 - 1.2 - 0.43 - clss70r_width
                    clss70r_react_per_lane = ll_react(1, clss70r_width, clss70r_min_pos, clss70r_max_pos, 0, Girder.no,
                                                      Girder.spacing, cantilever)
                    for k in range(Girder.no):
                        clss70r_react[k] += clss70r_react_per_lane[k]
                for k in range(Girder.no):
                    if clss70r_react[k] > react_max[k]:
                        react_max[k] = clss70r_react[k]

        return react_max


    def case2_clA(lcase_list):
        react_max = [0] * Girder.no
        for i in range(len(lcase_list)):
            clssA_lanes = lcase_list[i][0]
            clss70r_lanes = lcase_list[i][1]
            clssA_width = 1.8
            clssA_react = [0] * Girder.no
            if clssA_lanes == 0:
                pass
            else:
                clssA_min_pos = cb_left + 0.15 + 0.25
                if clss70r_lanes == 0:
                    clssA_max_pos = Deck.width - cb_right - 0.15 - 0.25 - 1.8 - (clssA_lanes - 1) * 3.5
                else:
                    clssA_max_pos = Deck.width - cb_right - 7.25 - (clss70r_lanes - 1) * 7 - 0.25 - 1.8 - \
                                    (clssA_lanes - 1) * 3.5
                clssA_react = ll_react(clssA_lanes, clssA_width, clssA_min_pos, clssA_max_pos, 3.5, Girder.no,
                                       Girder.spacing, cantilever)
                for k in range(Girder.no):
                    if clssA_react[k] > react_max[k]:
                        react_max[k] = clssA_react[k]
        return react_max


    def case2_70R(lcase_list):
        react_max = [0] * Girder.no
        for i in range(len(lcase_list)):
            clss70r_lanes = lcase_list[i][1]
            clss70r_width = 1.93
            clss70r_react = [0] * Girder.no
            if clss70r_lanes == 0:
                pass
            else:
                clss70r_react_per_lane = [0] * Girder.no
                for j in range(1, int(clss70r_lanes)):
                    clss70r_min_pos = Deck.width - cb_right - 7.25 - (j - 1) * 7 + 0.43 + 1.2
                    if j == 1:
                        clss70r_max_pos = Deck.width - cb_right - 1.2 - 0.43 - clss70r_width
                    elif j > 1:
                        clss70r_max_pos = Deck.width - cb_right - 7.25 - (j - 2) * 7 - 0.43 - clss70r_width
                    clss70r_react_per_lane = ll_react(1, clss70r_width, clss70r_min_pos, clss70r_max_pos, 0, Girder.no,
                                                      Girder.spacing, cantilever)
                    for k in range(Girder.no):
                        clss70r_react[k] += clss70r_react_per_lane[k]
                for k in range(Girder.no):
                    if clss70r_react[k] > react_max[k]:
                        react_max[k] = clss70r_react[k]
        return react_max


    a = case1_clA(lcase)
    b = case1_70R(lcase)
    c = case2_clA(lcase)
    d = case2_70R(lcase)


    class LiveLoad:
        outer_clA_case1 = max(a[0], a[Girder.no - 1])
        outer_70R_case1 = max(b[0], b[Girder.no - 1])
        outer_clA_case2 = max(c[0], c[Girder.no - 1])
        outer_70R_case2 = max(d[0], d[Girder.no - 1])
        inner_clA_case1 = 0
        for i in range(1, Girder.no - 2):
            if inner_clA_case1 < a[i]:
                inner_clA_case1 = a[i]
        inner_clA_case2 = 0
        for i in range(1, Girder.no - 2):
            if inner_clA_case2 < c[i]:
                inner_clA_case2 = c[i]
        inner_70r_case1 = 0
        for i in range(1, Girder.no - 2):
            if inner_70r_case1 < b[i]:
                inner_70r_case1 = b[i]
        inner_70r_case2 = 0
        for i in range(1, Girder.no - 2):
            if inner_70r_case2 < d[i]:
                inner_70r_case2 = d[i]
        sv_output = sv_load()
        sv_outer = max(sv_output[0], sv_output[Girder.no - 1])
        sv_inner = 0
        for i in range(1, Girder.no - 2):
            if sv_inner > sv_output[i]:
                sv_inner = sv_output[i]
        eff_span = Deck.length - Deck.long_cantilever * 2 - Girder.overhang * 2
        impact_factor_clA = 9 / (13.5 + eff_span) + 1
        congestion_factor_clA = 1
        i = 0
        while eff_span > codedata.span[i]:
            i += 1
        i -= 1
        congestion_factor_clA = (codedata.congestion_factor[i + 1] - codedata.congestion_factor[i]) / (
                    codedata.span[i + 1] - codedata.span[i]) * (eff_span - codedata.span[i]) + codedata.congestion_factor[i]
        load_factor_clA = impact_factor_clA * congestion_factor_clA
        impact_factor_70r = impact_factor_clA
        if eff_span < 23:
            impact_factor_70r = 1.25
        congestion_factor_70r = congestion_factor_clA
        load_factor_70r = impact_factor_70r * congestion_factor_70r


    # ######################## Load Summary
    def udl_bm(udl, location):
        overhang = Girder.overhang + Deck.long_cantilever
        bm = udl / 2 * (eff_span * location - overhang ** 2 - location ** 2)
        return bm


    def udl_sf(udl, location):
        sf = udl * (eff_span / 2 - location)
        return sf


    def ll_result(live_load, live_load_dist, location):
        overhang = Girder.overhang + Deck.long_cantilever
        loadlength = sum(live_load_dist)
        increment = 0.2
        totallength = eff_span + overhang * 2 + loadlength
        loadnumber = int((totallength - (totallength % increment)) / increment + 1)
        bm_max = 0
        sf_max = 0
        for i in range(loadnumber):
            bm = 0
            sf = 0
            loadoffset = increment * i - loadlength
            for j in range(len(live_load)):
                load = live_load[j]
                loadoffset += live_load_dist[j]
                if loadoffset < 0 or loadoffset > eff_span + 2 * overhang:
                    bm += 0
                    sf += 0
                elif loadoffset <= location + overhang:
                    bm += load / eff_span * (loadoffset - overhang) * (eff_span - location)
                    sf += load * (overhang + eff_span - loadoffset) / eff_span - load
                elif loadoffset > location + overhang:
                    bm += load * (overhang + eff_span - loadoffset) / eff_span * location
                    sf += load * (overhang + eff_span - loadoffset) / eff_span
            if bm > bm_max:
                bm_max = bm
            if abs(sf) > sf_max:
                sf_max = abs(sf)
        return [bm_max, sf_max]


    def bm_outer(location):
        girder_weight = udl_bm(Udl.girder_weight, location)
        deck_weight = udl_bm(Udl.deck_weight_outer, location)
        cb = udl_bm(Udl.sidl_outer, location)
        green_concrete = udl_bm(Udl.deck_weight_outer_green, location)
        shuttering = udl_bm(Udl.shutt_outer, location)
        const_ll = udl_bm(Udl.const_ll_outer, location)
        surfacing = udl_bm(Udl.surf_outer, location)
        wind_const = udl_bm(Udl.wind_const_outer, location)
        wind_no_ll = udl_bm(Udl.wind_noll_outer, location)
        wind_ll = udl_bm(Udl.wind_ll_outer, location)
        seis_no_ll = udl_bm(Udl.seis_noll_outer, location)
        live_load_case1 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[0] * LiveLoad.load_factor_clA * LiveLoad.outer_clA_case1 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[0] * LiveLoad.load_factor_70r * LiveLoad.outer_70R_case1 / 100
        live_load_case2 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[0] * LiveLoad.load_factor_clA * LiveLoad.outer_clA_case2 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[0] * LiveLoad.load_factor_70r * LiveLoad.outer_70R_case2 / 100
        live_load = max(live_load_case1, live_load_case2)
        live_load_sv = ll_result(codedata.class_sv, codedata.class_sv_dist, location)[0] * LiveLoad.sv_outer / 100
        live_load_seis = live_load * 0.2 * Loads.seis_coeff
        load = [girder_weight, deck_weight, cb, green_concrete, shuttering, const_ll, surfacing, wind_const, wind_no_ll,
                wind_ll, seis_no_ll, live_load, live_load_sv, live_load_seis]
        return load


    def bm_inner(location):
        girder_weight = udl_bm(Udl.girder_weight, location)
        deck_weight = udl_bm(Udl.deck_weight_inner, location)
        cb = udl_bm(Udl.sidl_inner, location)
        green_concrete = udl_bm(Udl.deck_weight_inner_green, location)
        shuttering = udl_bm(Udl.shutt_inner, location)
        const_ll = udl_bm(Udl.const_ll_inner, location)
        surfacing = udl_bm(Udl.surf_inner, location)
        wind_const = udl_bm(Udl.wind_const_inner, location)
        wind_no_ll = udl_bm(Udl.wind_noll_inner, location)
        wind_ll = udl_bm(Udl.wind_ll_inner, location)
        seis_no_ll = udl_bm(Udl.seis_noll_inner, location)
        live_load_case1 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[0] * LiveLoad.load_factor_clA * LiveLoad.inner_clA_case1 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[0] * LiveLoad.load_factor_70r * LiveLoad.inner_70r_case1 / 100
        live_load_case2 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[0] * LiveLoad.load_factor_clA * LiveLoad.inner_clA_case2 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[0] * LiveLoad.load_factor_70r * LiveLoad.inner_70r_case2 / 100
        live_load = max(live_load_case1, live_load_case2)
        live_load_sv = ll_result(codedata.class_sv, codedata.class_sv_dist, location)[0] * LiveLoad.sv_inner / 100
        live_load_seis = live_load * 0.2 * Loads.seis_coeff
        load = [girder_weight, deck_weight, cb, green_concrete, shuttering, const_ll, surfacing, wind_const, wind_no_ll,
                wind_ll, seis_no_ll, live_load, live_load_sv, live_load_seis]
        return load


    def sf_outer(location):
        girder_weight = udl_sf(Udl.girder_weight, location)
        deck_weight = udl_sf(Udl.deck_weight_outer, location)
        cb = udl_sf(Udl.sidl_outer, location)
        green_concrete = udl_sf(Udl.deck_weight_outer_green, location)
        shuttering = udl_sf(Udl.shutt_outer, location)
        const_ll = udl_sf(Udl.const_ll_outer, location)
        surfacing = udl_sf(Udl.surf_outer, location)
        wind_const = udl_sf(Udl.wind_const_outer, location)
        wind_no_ll = udl_sf(Udl.wind_noll_outer, location)
        wind_ll = udl_sf(Udl.wind_ll_outer, location)
        seis_no_ll = udl_sf(Udl.seis_noll_outer, location)
        live_load_case1 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[1] * LiveLoad.load_factor_clA * LiveLoad.outer_clA_case1 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[1] * LiveLoad.load_factor_70r * LiveLoad.outer_70R_case1 / 100
        live_load_case2 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[1] * LiveLoad.load_factor_clA * LiveLoad.outer_clA_case2 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[1] * LiveLoad.load_factor_70r * LiveLoad.outer_70R_case2 / 100
        live_load = max(live_load_case1, live_load_case2)
        live_load_sv = ll_result(codedata.class_sv, codedata.class_sv_dist, location)[1] * LiveLoad.sv_outer / 100
        live_load_seis = live_load * 0.2 * Loads.seis_coeff
        load = [girder_weight, deck_weight, cb, green_concrete, shuttering, const_ll, surfacing, wind_const, wind_no_ll,
                wind_ll, seis_no_ll, live_load, live_load_sv, live_load_seis]
        return load


    def sf_inner(location):
        girder_weight = udl_sf(Udl.girder_weight, location)
        deck_weight = udl_sf(Udl.deck_weight_inner, location)
        cb = udl_sf(Udl.sidl_inner, location)
        green_concrete = udl_sf(Udl.deck_weight_inner_green, location)
        shuttering = udl_sf(Udl.shutt_inner, location)
        const_ll = udl_sf(Udl.const_ll_inner, location)
        surfacing = udl_sf(Udl.surf_inner, location)
        wind_const = udl_sf(Udl.wind_const_inner, location)
        wind_no_ll = udl_sf(Udl.wind_noll_inner, location)
        wind_ll = udl_sf(Udl.wind_ll_inner, location)
        seis_no_ll = udl_sf(Udl.seis_noll_inner, location)
        live_load_case1 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[1] * LiveLoad.load_factor_clA * LiveLoad.inner_clA_case1 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[1] * LiveLoad.load_factor_70r * LiveLoad.inner_70r_case1 / 100
        live_load_case2 = ll_result(codedata.classA, codedata.classA_dist,
                                    location)[1] * LiveLoad.load_factor_clA * LiveLoad.inner_clA_case2 / 100 + \
                          ll_result(
                              codedata.class_70r, codedata.class_70r_dist,
                              location)[1] * LiveLoad.load_factor_70r * LiveLoad.inner_70r_case2 / 100
        live_load = max(live_load_case1, live_load_case2)
        live_load_sv = ll_result(codedata.class_sv, codedata.class_sv_dist, location)[1] * LiveLoad.sv_inner / 100
        live_load_seis = live_load * 0.2 * Loads.seis_coeff
        load = [girder_weight, deck_weight, cb, green_concrete, shuttering, const_ll, surfacing, wind_const, wind_no_ll,
                wind_ll, seis_no_ll, live_load, live_load_sv, live_load_seis]
        return load


    def factored_load(load, factor_list):
        max_load = 0
        for i in factor_list:
            p = 0
            for j in range(len(load)):
                p += load[j] * codedata.factor[i][j]
            if p > max_load:
                max_load = p
        return max_load


    def bm_stage1(load):
        factor_list = [0]
        factored_load_value = factored_load(load, factor_list)
        bending_strength = PlateGirder.bending_strength_stage1
        ur = factored_load_value / bending_strength
        return ur


    def sf_stage1(load):
        factor_list = [0]
        factored_load_value = factored_load(load, factor_list)
        shear_strength = PlateGirder.shear_strength
        ur = factored_load_value / shear_strength
        return ur


    def bm_stage2(load):
        factor_list = [1]
        factored_load_value = factored_load(load, factor_list)
        bending_strength = PlateGirder.bending_strength_stage2
        ur = factored_load_value / bending_strength
        return ur


    def sf_stage2(load):
        factor_list = [1]
        factored_load_value = factored_load(load, factor_list)
        shear_strength = PlateGirder.shear_strength
        ur = factored_load_value / shear_strength
        return ur


    def bm_outer_stage3(load):
        factor_list = [1]
        factored_load_value = factored_load(load, factor_list)
        bending_strength = CompGirder.outer_bending_strength
        ur = factored_load_value / bending_strength
        return ur


    def bm_inner_stage3(load):
        factor_list = [1]
        factored_load_value = factored_load(load, factor_list)
        bending_strength = CompGirder.inner_bending_strength
        ur = factored_load_value / bending_strength
        return ur


    def sf_stage3(load):
        factor_list = [1]
        factored_load_value = factored_load(load, factor_list)
        shear_strength = PlateGirder.shear_strength
        ur = factored_load_value / shear_strength
        return ur


    ur_outer_bm_stage1 = []
    ur_outer_bm_stage2 = []
    ur_outer_bm_stage3 = []

    ur_outer_sf_stage1 = []
    ur_outer_sf_stage2 = []
    ur_outer_sf_stage3 = []

    ur_inner_bm_stage1 = []
    ur_inner_bm_stage2 = []
    ur_inner_bm_stage3 = []

    ur_inner_sf_stage1 = []
    ur_inner_sf_stage2 = []
    ur_inner_sf_stage3 = []

    divisions = 4
    segment = eff_span / divisions / 2
    for i in range(divisions + 1):
        ur_outer_bm_stage1.append(bm_stage1(bm_outer(segment * i)))
        ur_outer_bm_stage2.append(bm_stage2(bm_outer(segment * i)))
        ur_outer_bm_stage3.append(bm_outer_stage3(bm_outer(segment * i)))

        ur_outer_sf_stage1.append(sf_stage1(sf_outer(segment * i)))
        ur_outer_sf_stage2.append(sf_stage2(sf_outer(segment * i)))
        ur_outer_sf_stage3.append(sf_stage3(sf_outer(segment * i)))

        ur_inner_bm_stage1.append(bm_stage1(bm_inner(segment * i)))
        ur_inner_bm_stage2.append(bm_stage2(bm_inner(segment * i)))
        ur_inner_bm_stage3.append(bm_inner_stage3(bm_inner(segment * i)))

        ur_inner_sf_stage1.append(sf_stage1(sf_inner(segment * i)))
        ur_inner_sf_stage2.append(sf_stage2(sf_inner(segment * i)))
        ur_inner_sf_stage3.append(sf_stage3(sf_inner(segment * i)))
    ur_list = [ur_outer_bm_stage1,
               ur_outer_bm_stage2,
               ur_outer_bm_stage3,
               ur_outer_sf_stage1,
               ur_outer_sf_stage2,
               ur_outer_sf_stage3,
               ur_inner_bm_stage1,
               ur_inner_bm_stage2,
               ur_inner_bm_stage3,
               ur_inner_sf_stage1,
               ur_inner_sf_stage2,
               ur_inner_sf_stage3]
    desired_ur_list = [0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95]
    status = 0
    for i in range(12):
        if max(ur_list[i]) > desired_ur_list[i]:
            status += 1
    return status


ffg = [0.5, 0.75, 1.0]
for j in range(100):
    for i in ffg:
        print(design_run(ffg))