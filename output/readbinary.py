import numpy as np

dens_binary = np.fromfile("binary_dens_G512_BS2500_FAC1_S1000000.dat",  dtype=np.float32)

print(len(dens_binary))
#8519633
N = 8388608
print(dens_binary[N:N+10])

# i, dens = np.loadtxt("txt_dens_G512_BS2500_FAC1_S1000000.dat", usecols=(0,1),unpack=True)
# it, denst = np.loadtxt("binary_txt_dens_G512_BS2500_FAC1_S1000000.dat", usecols=(0,1),unpack=True)

# print(dens[0:10])
# print(denst[0:10])


# print(dens_binary[-10:])
# print(dens[-10:])
# print(denst[-10:])