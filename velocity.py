import numpy as np


def boundary(domain,bc,m,n):
    #Vertical Wall
    for j in range(n):
        domain[j][0] = bc[0]        #Leftside Wall
        domain[j][m-1] = bc[1]      #Rightside Wall
        
    #Horizontal Wall
    for i in range(m):
        domain[0][i] = bc[2]        #Bottomside Wall
        domain[n-1][i] = bc[3]      #Topside Wall

    return domain

def u_velocity(psi,u,dy,m,n):    
    for j in range(1,n-1):
        for i in range(1,m-1):
            u[j][i] = (psi[j][i+1]-psi[j][i-1])/(2*dy)
    return u

def v_velocity(psi,v,dx,m,n):    
    for j in range(1,n-1):
        for i in range(1,m-1):
            v[j][i] = (psi[j+1][i]-psi[j-1][i])/(2*dx)
    return v