import numpy as np

def boundary(omega,psi,u,v,dx,dy,m,n):
    for j in range(0,n):
        omega[j][0] = (-2*(psi[j][1])/(dx*dx))-(2*u[j][0]/dx)           #Left wall
        omega[j][m-1] = (-2*(psi[j][n-1])/(dx*dx))-(2*u[j][n-1]/dx)     #Right wall

    for i in range(0,m):
        omega[0][i] = (-2*(psi[1][i])/(dy*dy))-(2*u[0][i]/dx)           #Top wall
        omega[n-1][i] = (-2*(psi[m-1][i])/(dy*dy))-(2*u[m-1][i]/dx)     #Bottom wall

    return omega


def transposrt(n_domain,o_domain,psi,dx,dy,dt,m,n,nu):
    for i in range(1,m-1):
        for j in range(1,n-1):
            n_domain[j][i] = o_domain[j][i] +\
                dt*(((psi[j][i+1]-psi[j][i-1])*(o_domain[j+1][i]-o_domain[j-1][i])/(4*dx*dx))-\
                    ((psi[j+1][i]-psi[j-1][i])*(o_domain[j][i+1]-o_domain[j][i-1])/(4*dy*dy))+\
                        (nu*(o_domain[j][i+1]+o_domain[j][i-1]+o_domain[j+1][i]+o_domain[j-1][i]-4*o_domain[j][i])/(dx*dx)))
    return n_domain