####################################################################
#                           F3OF
####################################################################

import numpy as np
import matplotlib.pyplot as plt
from inwavepp.mesh.LoadMesh import loadVTK
from inwavepp.mesh.ProcessMesh import invertNormals
from innosea_nemoh_library.mesh.importMesh import loadMSH, writeL12
from innosea_nemoh_library.mesh.importMesh import importMSH
from innosea_nemoh_library.mesh.customMesh import customMesh
from innosea_nemoh_library.nemohDriver.Nemoh import HydroSolver
from innosea_nemoh_library.motion.mechanicalSolver import mechanicalSolver
from innosea_nemoh_library.plot.plotResults import PostProcessing
from innosea_nemoh_library.loadData.loadNemoh import loadnBodies, loadnemohcal, loadw, loadDir, loadFe, loadMa, loadMa2, loadB,\
    loadnNodesPanels
from inwavepp.utils.PostProcessingParameters import MARKERS_REF, MARKERS_INW, xlabel_freq
from inwavepp.plotResults.PlotResults import getnewfigure, showplots, plotXY, saveplots

# Convert .msh files to .dat format

# base_vtk = "Mesh\\base_eq.vtk"
# base_xyz, base_NN  = loadVTK(base_vtk)
# base_dat = "Mesh\\base.dat"
# writeL12(base_xyz, base_NN, base_dat)
#
# flap_1_vtk = "Mesh\\flap 1_eq.vtk"
# flap_1_xyz, flap_1_NN = loadVTK(flap_1_vtk)
# flap_1_dat = "Mesh\\flap_1.dat"
# writeL12(flap_1_xyz, flap_1_NN, flap_1_dat)
#
# flap_2_vtk = "Mesh\\flap 2_eq.vtk"
# flap_2_xyz, flap_2_NN = loadVTK(flap_2_vtk)
# flap_2_dat = "Mesh\\flap_2.dat"
# writeL12(flap_2_xyz, flap_2_NN, flap_2_dat)

# Nemoh
#HydroSolver()

nBodies = loadnBodies()
totalndofs, nw, w, ndir, Dir, DirKochin, FSEParam, PressureParam = loadnemohcal()
#modulesFe, phasesFe = loadFe()
Ma = loadMa()
Bw = loadB()

# convert to inw frame
def skew(x):
    return np.array([[0, -x[2], x[1]],
                     [x[2], 0, x[0]],
                     [x[1], x[0], 0]])

R = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])

Pb = np.array([0, 0, -9])
P1 = np.array([-12.5, 0, -4.25])
P2 = np.array([12.5, 0, -4.25])

Pbskew = skew(Pb)
P1skew = skew(P1)
P2skew = skew(P2)

Xbfor = np.array([[R, R * Pbskew.transpose()], [np.zeros((3, 3)), R]])
Xbback = np.array([[R, np.zeros((3, 3))], [Pbskew * R, R]])

X1for = np.array([[R, R * P1skew.transpose()], [np.zeros((3, 3)), R]])
X1back = np.array([[R, np.zeros((3, 3))], [P1skew * R, R]])

X2for = np.array([[R, R * P2skew.transpose()], [np.zeros((3, 3)), R]])
X2back = np.array([[R, np.zeros((3, 3))], [P2skew * R, R]])

#print(Xbfor)
#print(Xbfor)
#print(X1for)
#print(X1for)
#print(X2for)
#print(X2for)

Ma_11 = Ma[0, 0, :] + Ma[0, 6, :] + Ma[0, 12, :] + Ma[6, 0, :] + Ma[6, 6, :]\
        + Ma[6, 12, :] + Ma[12, 0, :] + Ma[12, 6, :] + Ma[12, 12, :]
Ma_33 = Ma[2, 2, :] + Ma[2, 8, :] + Ma[2, 14, :] + Ma[8, 2, :] + Ma[8, 8, :]\
        + Ma[8, 14, :] + Ma[14, 2, :] + Ma[14, 8, :] + Ma[14, 14, :]
Ma_55 = Ma[4, 4, :] + Ma[4, 10, :] + Ma[4, 16, :] + Ma[10, 4, :] + Ma[10, 10, :]\
        + Ma[10, 16, :] + Ma[16, 4, :] + Ma[16, 10, :] + Ma[16, 16, :]
B_11 = Bw[0, 0, :] + Bw[0, 6, :] + Bw[0, 12, :] + Bw[6, 0, :] + Bw[6, 6, :] \
        + Bw[6, 12, :] + Bw[12, 0, :] + Bw[12, 6, :] + Bw[12, 12, :]
B_33 = Bw[2, 2, :] + Bw[2, 8, :] + Bw[2, 14, :] + Bw[8, 2, :] + Bw[8, 8, :] \
        + Bw[8, 14, :] + Bw[14, 2, :] + Bw[14, 8, :] + Bw[14, 14, :]
B_55 = Bw[4, 4, :] + Bw[4, 10, :] + Bw[4, 16, :] + Bw[10, 4, :] + Bw[10, 10, :] \
        + Bw[10, 16, :] + Bw[16, 4, :] + Bw[16, 10, :] + Bw[16, 16, :]

Ma_pitch00 = Ma[4, 4, :]
Ma_pitch11 = Ma[4, 10, :]
Ma_pitch22 = Ma[16, 16, :]

Bw_pitch00 = Bw[4, 4, :]
Bw_pitch11 = Bw[10, 10, :]
Bw_pitch22 = Bw[16, 16, :]

plt.plot(w, B_33, label='NEMOH')
plt.show()

np.save('prints/B_data', Bw)
np.save('prints/Ma_data', Ma)
#handle = getnewfigure()
#plotXY(w, Ma_11, handle=handle, marker=MARKERS_INW[0], yname='NEMOH', xlabel=xlabel_freq,
#       ylabel="RAO", title= 'Ma_11', style='sci')
#handle = getnewfigure()
#plotXY(w, Ma_33, handle=handle, marker=MARKERS_INW[0], yname='NEMOH', xlabel=xlabel_freq,
#       ylabel="RAO", title='Ma_33', style='sci')
#handle = getnewfigure()
#plotXY(w, Ma_55, handle=handle, marker=MARKERS_INW[0], yname='NEMOH', xlabel=xlabel_freq,
#       ylabel="RAO", title='Ma_55', style='sci')
#handle = getnewfigure()
#plotXY(w, B_11, handle=handle, marker=MARKERS_INW[0], yname='NEMOH', xlabel=xlabel_freq,
#       ylabel="RAO", title='B_11', style='sci')
#handle = getnewfigure()
#plotXY(w, B_33, handle=handle, marker=MARKERS_INW[0], yname='NEMOH', xlabel=xlabel_freq,
#       ylabel="RAO", title='B_33', style='sci')
#handle = getnewfigure()
#plotXY(w, B_55, handle=handle, marker=MARKERS_INW[0], yname='NEMOH', xlabel=xlabel_freq,
#       ylabel="RAO", title='B_55', style='sci')

#
# plt.plot(w, Ma_11, label='NEMOH')
# #plt.show()
# plt.savefig("Ma_11.png")
# plt.plot(w, Ma_33, label='NEMOH')
# #plt.show()
# plt.savefig("Ma_33.png")
#plt.plot(w, Ma_55, label='NEMOH')
#plt.show()
# plt.savefig("Ma_55.png")
#

#
# plt.plot(w, B_11, label='NEMOH')
# #plt.show()
# plt.savefig("B_11.png")
# plt.plot(w, B_33, label='NEMOH')
# #plt.show()
# plt.savefig("B_33.png")

# plt.savefig("B_55.png")

#saveplots('figures/')

#nNodes, nPanels = loadnNodesPanels("\\flap_2")

#print "nNodes = ", nNodes
#print " nPanels = ", nPanels

# Mechanical solver
#mechanicalSolver()

# Plots
#PostProcessing(save = True,show = True)



# # plt.plot(RAOs_T, RAOs_HOTINT, label='HOTINT (FeFreq)', linestyle='--', marker='o', color='r')
# # plt.plot(RAOs_T, RAOs_MARIN, label='MARIN', linestyle='--', marker='s', color='k')
# # plt.plot(RAOs_T, RAOs_INWAVE_dw, label='InWave (0.02<w<6.0, dw=0.02, 40s IRF)', linestyle='--', marker='^',
# color='g')
# # plt.plot(RAOs_T, RAOs_INWAVE_IRF, label='InWave (0.05<w<12.0, dw=0.05, 10s IRF)', linestyle='--', marker='h',
# # color='b')
# # plt.xlim(0, 12)
# # plt.xlabel('Tp (s)')
# # plt.ylabel('Heave RAO')
# # #plt.title('title')
# # plt.grid(True)
# # plt.legend()
# # #plt.savefig("test.png")