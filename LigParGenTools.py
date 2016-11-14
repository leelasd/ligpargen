from math import *
from random import randint, random
import numpy as np
############################################


def trans(coos, qi):
    rcoos = np.zeros((len(coos), len(qi)))
    for i in range(len(coos)):
        rcoos[i][0] = coos[i][0] + qi[0]
        rcoos[i][1] = coos[i][1] + qi[1]
        rcoos[i][2] = coos[i][2] + qi[2]
    return(rcoos)
#################### CENTER OF COOS #############################


def CenOfCoos(atoms, coos):
    com = np.zeros(3)
    for i in range(0, len(atoms)):
        com[0] += coos[i][0]
        com[1] += coos[i][1]
        com[2] += coos[i][2]
    ncom = [i / len(atoms) for i in com]
    rad = []
    for i in range(0, len(atoms)):
        coos[i][0] = coos[i][0] - ncom[0]
        coos[i][1] = coos[i][1] - ncom[1]
        coos[i][2] = coos[i][2] - ncom[2]
        rad.append(sqrt(coos[i][0]**2 + coos[i][1]**2 + coos[i][2]**2))
    return(coos, ncom, max(rad))
####################MOMENTS OF INERTIA #############################


def rotzyx(coos):
    phi = random() * 2 * np.pi
    csphi = np.cos(phi)
    snphi = np.sin(phi)
    chi = random() * 2 * np.pi
    cschi = np.cos(chi)
    snchi = np.sin(chi)
    cstta = (random() * 2) - 1
    sntta = np.sin(np.arccos(cstta))
    RXX = (cstta * csphi * cschi) - (snphi * snchi)
    RXY = -cstta * csphi * snchi - snphi * cschi
    RXZ = sntta * csphi
    RYX = cstta * snphi * cschi + csphi * snchi
    RYY = (-cstta * snphi * snchi) + (csphi * cschi)
    RYZ = sntta * snphi
    RZX = -sntta * cschi
    RZY = sntta * snchi
    RZZ = cstta
    rcoos = np.zeros((len(coos), 3))
    for i in range(len(coos)):
        rcoos[i][0] = coos[i][0] * RXX + coos[i][1] * RXY + coos[i][2] * RXZ
        rcoos[i][1] = coos[i][0] * RYX + coos[i][1] * RYY + coos[i][2] * RYZ
        rcoos[i][2] = coos[i][0] * RZX + coos[i][1] * RZY + coos[i][2] * RZZ
    return (rcoos)
########################## FOR READING BOXMAKER COOS ######################


def read_pdb(fname):
    fcon = open(fname).readlines()
    atoms = []
    coos = []
    pdb_lines = {}
    for line in fcon:
        if ('ATOM' in line) or ('HETATM' in line):
            atoms.append(line[12:16].strip())
            coos.append(list(map(float, line[32:54].split())))
            # pdb_lines[line[12:16].strip()]=list(map(float,line[32:54].split()))
    return atoms, coos

########################## FOR READING BOXMAKER COOS ######################


def pdb_lines(atoms, coos, ID, resid='UNK'):
    lines = []
    num = (ID - 1) * len(atoms)
    for (i, j) in zip(atoms, coos):
        num += 1
        lines.append('%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f\n' %
                     ('ATOM', num, i, resid, ID, j[0], j[1], j[2]))
    lines.append('TER \n')
    return lines

########################## FOR READING BOXMAKER COOS ######################


def BOX_MAKER(pdb_file, BOX_SIZE):
    atMOL, csMOL = read_pdb(pdb_file)
    csMOL, comMOL, radMOL = CenOfCoos(atMOL, csMOL)
    maxX = 1.5 + (max([c[0] for c in csMOL]) - min([c[0] for c in csMOL]))
    maxY = 1.5 + (max([c[1] for c in csMOL]) - min([c[1] for c in csMOL]))
    maxZ = 1.5 + (max([c[2] for c in csMOL]) - min([c[2] for c in csMOL]))
    (nx, ny, nz) = (int(BOX_SIZE / maxX),
                    int(BOX_SIZE / maxY), int(BOX_SIZE / maxZ))
    gx = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=nx)
    gy = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=ny)
    gz = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=nz)

    BOX = {}
    ID = 0
    total_lines = []
    for XC in gx:
        for YC in gy:
            for ZC in gz:
                csMOL = rotzyx(csMOL)
                csF = trans(csMOL, [XC, YC, ZC])
                BOX[ID] = {'ATS': atMOL, 'CS': csF}
                total_lines += pdb_lines(BOX[ID]['ATS'], BOX[ID]['CS'], ID + 1)
                ID = ID + 1

    total_mols = ID
    ofile = open('box_%d_' % (int(BOX_SIZE)) + pdb_file, 'w+')
    for l in total_lines:
        ofile.write('%s' % l)
    ofile.close()
    return None
########################## FOR READING BOXMAKER COOS ######################
