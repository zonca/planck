import pointingtools
from planck import Planck
pl=Planck.Planck()
s=pointingtools.SiamAngles(False)
lfi=pl.inst['LFI']
lines=[]
lines.append('2011-05-23 TestEnv JupE 1.0\n')
for ch in lfi.ch:
    lines.append(ch.tag + '\n')
    ss=s.get(ch).T
    for row in ss:
        lines.append('%.16e %.16e %.16e\n' % tuple(row))
a = open('siam_instrument_jupE1.0.txt', 'w')
a.writelines(lines)
a.close()
