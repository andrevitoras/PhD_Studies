"""
Sharma, V., Khanna, S., Nayak, J.K., Kedare, S.B., 2016.
Effects of shading and blocking in compact linear fresnel reflector field.
Energy 94, 633–653. https://doi.org/10.1016/J.ENERGY.2015.10.098
"""
from numpy import array

if __name__ == '__main__':

    def get_sharma_optimum_data(p_to_w: float,
                                n=25, Hr_to_nw=0.7):
        """
        This function is based on the optimum data reported in Table 3, p. 131.

        For this case, as detailed in section 4.2,
        the authors state that a fixed value of n = 25 is considered, where n is the number of primaries.

        This optimum calculations are based on their argument that for n > 8,
        "any combinations of H/w and n corresponding to a fixed H/nw yields the same values of annual factors".
        These calculations were carried out by considering p/w = 1.42 and \Omega = 0.

        :param p_to_w: mirrors pitch to width ratio.
        :param Hr_to_nw: receiver height to mirror width times number of mirrors. Hr/nw = 0.7, as defined by the authors
        :param n: number of primary mirrors. n = 25, as defined by the authors..

        :return: the filling factor and the receiver height to primary width ratio
        """

        # It is easy to derive the following expressions
        FF = n / ((n-1)*p_to_w + 1)  # the filling factor

        Hr_to_wp = Hr_to_nw * FF  # receiver height to primary width ratio

        return FF, Hr_to_wp

    optimum_p_to_w = [1.65, 1.22, 1.75, 1.39, 1.86, 1.63]

    optimum_data = array([get_sharma_optimum_data(x) for x in optimum_p_to_w])

    optimum_FF = optimum_data.T[0].round(3)
    optimum_Hr_to_wp = optimum_data.T[1].round(3)

    print(f'The optimum values of Filling Factor are {optimum_FF}!\n')
    print(f'The optimum values of Hr/Wp are {optimum_Hr_to_wp}!\n')

