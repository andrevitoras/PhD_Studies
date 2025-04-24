

if __name__ == '__main__':

    # Data from Table 2, p.195 -- base values before optimization
    number_primaries = 25.0
    mirror_width = 0.6
    mirrors_orientation = 'NS'
    primary_field_width = 21.0

    s = (primary_field_width - mirror_width) / (number_primaries - 1)
    print(f'Before optimization, the shift between primary mirrors is {s}!')

    # Calculations ####################################################################################################

    # Values that authors consider as the optimum one
    # they are based on the useful energy efficiency metric
    # p. 199 - right side column
    FF = 0.72  # filling factor
    Hr_to_wp = 0.5

    # This work consider uniform primary fields.
    # the number of mirrors and mirror width are constant, so that varying the filling factor changes the distance
    # between mirrors and the corresponding primary field width

    wp = (number_primaries * mirror_width) / FF
    s = (wp - mirror_width) / (number_primaries - 1)

    s_to_w = s / mirror_width
    print(f'For the optimum geometry reported by the authors, the s/w ratio is {round(s/mirror_width, 3)}!')
