import math

#
#   Useful regular expression for QE
#
decimal=r'-?\d+\.\d*|-?\.\d+'
real_sn=r'[+\-]?\d*\.\d*[E][+\-]\d\d?'
r_or_d =r'('+real_sn + '|' + decimal + ')'


def bose(E,T):
    #
    # Bose function
    if T < 1e-10:
        return 0.0
    else:
        return 1.0/(math.exp(E/T)-1.0)
