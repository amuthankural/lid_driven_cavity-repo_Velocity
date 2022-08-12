import sys
from os import stat, stat_result
from datetime import datetime, time
import json
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import shape

from numpy.lib.index_tricks import c_

np.set_printoptions(precision=2)
np.set_printoptions(threshold=sys.maxsize)

def prog_run():
    with open('./Input.json','r') as f:
        input           = json.load(f)

    """-----------Imported data----------"""
    l_i                 = input["lx_i"]                 #Length of the domain
    b_j                 = input["by_j"]                 #Breadth of the domain
    n_i                 = input["n_i"]                  #Number of grid points along 'x'
    m_j                 = input["m_j"]                  #Number of grid points along 'y'
    Re                  = input["Re"]                   #Reynolds Number
    bc_u                = [input["u_L"],input["u_R"],input["u_B"],input["u_T"]]
    bc_v                = [input["v_L"],input["v_R"],input["v_B"],input["v_T"]]
    criteria            = input["criteria"]
    status              = False



    """-----------Calculated data----------"""
    dx                  = l_i/(n_i-1)                   #Grid length along 'x' direction
    dy                  = b_j/(m_j-1)                   #Grid length along 'y' direction
    beta                = dx/dy
    dim                     = (n_i,m_j)
    
    """---------------Domains--------------"""
    om_old                  = discretize(n_i,m_j)
    om_domain               = discretize(n_i,m_j)
    si_domain               = discretize(n_i,m_j)
    u_domain                = discretize(n_i,m_j)
    v_domain                = discretize(n_i,m_j)
    #vel_domain              = np.zeros([n_i, m_j, 2], dtype= float)
    #print(vel_domain.shape)
    print((dx,dy))

    u_domain = vel_bound_condition(u_domain,bc_u,dim)
    v_domain = vel_bound_condition(v_domain,bc_v,dim)
    om_domain = om_calc(om_domain,u_domain,v_domain,dx,dy,dim)
    si_domain = psi_calc(om_domain,si_domain,dx,beta,dim)
    u_domain = u_calc(si_domain,u_domain,dy,dim)
    v_domain = v_calc(si_domain,v_domain,dx,dim)
    print(u_domain)
    


"""
    print("U - Velocity:\n",u_domain)
    print("\nV - Velocity:\n",v_domain)

    while (status == False):
        iteration = 1
        om_calc(om_domain,u_domain,v_domain,dim)
        psi_calc(om_domain,si_domain,beta)
        u_calc(si_domain,u_domain,dy)
        v_calc(si_domain,v_domain,dx)
        np.copyto(om_old,om_domain)
        status = converge(om_old,om_domain,criteria)
        iteration += 1
"""

def discretize(n_i,m_j):
    return(np.zeros([n_i, m_j], dtype= float))

def vel_bound_condition(domain,bc,dim):
    for i in range(dim[0]):
        domain[i][0] = bc[0]
        domain[i][dim[1]-1] = bc[1]
        
    for j in range(dim[1]):
        domain[0][j] = bc[2]
        domain[dim[0]-1][j] = bc[3]
        
    return domain

#Central difference scheme
def om_calc(omega,u,v,dx,dy,dim):
    print("dimension: ",dim,"\n")
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            omega[i][j] = ((v[i][j+1] - 2*v[i][j] + v[i][j-1])/(dx*dx)) - ((u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dy*dy))
            #print("iteration: ",i,j,"\n")
    return omega


def vorticity(om,u,v,dx,dy,nu,dim):
    
    return om


#Elliptic equation (Poisson equation)
def psi_calc(om,si,dx,beta,dim):
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            si[i][j] = ((si[i][j+1]+si[i][j-1]+\
                ((beta*beta)*(si[i+1][j]+si[i-1][j]))+\
                    (om[i][j]*dx*dx))/\
                        (2*(1+beta*beta)))
    return si

def u_calc(si,u,dy,dim):    
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            u[i][j] = (si[i+1][j]-si[i-1][j])/(2*dy)
    return u

def v_calc(si,v,dx,dim):    
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            v[i][j] = (si[i][j+1]-si[i][j-1])/(2*dx)
    return v



def converge(om_old,om_new,criteria):
    convergence = np.subtract(om_old,om_new)
    if (np.amax(convergence,axis = None) <= criteria):
        return True
    else:
        return False

prog_run()
