import numpy as np

def poisson(omega,psi,dx,m,n):
    for j in range(1,n-1):
        for i in range(1,m-1):
            psi[j][i] = ((dx*dx*omega[j][i])+psi[j][i+1]+psi[j][i-1]+psi[j+1][i]+psi[j-1][i])/4
    return psi

