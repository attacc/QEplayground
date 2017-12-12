from QEplayground.units import *


mass=12.0107*amu2au

#freq=33.126932*thz2cm1*cm1toeV*ev2ha
freq=0.137*0.037

der=992.12*(0.529177249)**2

h=1.0

const=h/(mass*freq)

print(const/8)
print(const*ha2ev/2)
#print(const*der*ha2ev)
