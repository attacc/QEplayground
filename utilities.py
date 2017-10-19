#
#   Useful regular expression for QE
#
decimal=r'-?\d+\.\d*|-?\.\d+'
real_sn=r'[+\-]?\d*\.\d*[E][+\-]\d\d?'
r_or_d =r'('+real_sn + '|' + decimal + ')'

#
#   Units
# 
AngToAu=1.8897261
