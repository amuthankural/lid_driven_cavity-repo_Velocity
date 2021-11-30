import sys
import logging
from os import stat, stat_result
from datetime import datetime, time
import json
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import shape

from numpy.lib.index_tricks import c_

np.set_printoptions(precision=2)
np.set_printoptions(threshold=sys.maxsize)



logging.basicConfig(filename='ftcs.log', encoding='utf-8', level=logging.DEBUG,\
    format = '%(asctime)s:%(name)s:%(message)s')

def prog_run():
    with open('./input.json','r') as f:
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
    nu                  = input["nu"]
    dt                  = input["dt"]
    status              = False
    logging.debug('Data imported from JSON file')


    """-----------Calculated data----------"""
    dx                  = l_i/(n_i-1)                   #Grid length along 'x' direction
    dy                  = b_j/(m_j-1)                   #Grid length along 'y' direction
    beta                = dx/dy
    dim                 = (n_i,m_j)
    
    """---------------Domains--------------"""
    n_om_domain             = discretize(n_i,m_j)
    om_domain               = discretize(n_i,m_j)
    si_domain               = discretize(n_i,m_j)
    u_domain                = discretize(n_i,m_j)
    v_domain                = discretize(n_i,m_j)
    logging.debug('Domain arrays created')
    
    u_domain = vel_bound_condition(u_domain,bc_u,dim)
    v_domain = vel_bound_condition(v_domain,bc_v,dim)
    om_domain = om_calc(om_domain,u_domain,v_domain,dx,dy,dim,criteria)
    si_domain = psi_calc(om_domain,si_domain,dx,beta,dim)


    iteration = 1
    logging.debug('Stepping into iterations')
    while (status == False):
        logging.debug('-------------------------------------------------------------------------------------------------')
        logging.debug('Iteration step: {}'.format(iteration))
        logging.debug('U domain at iteration start: \n{}\n'.format(u_domain))
        logging.debug('V domain at iteration start: \n{}\n'.format(v_domain))
        logging.debug('Initial Vorticity domain: \n{}\n'.format(om_domain))
        logging.debug('Initial StreamFunction domain: \n{}\n'.format(si_domain))
        om_domain = vorticity_bound(om_domain,si_domain,u_domain,v_domain,dx,dy,dim)
        logging.debug('Vorticity domain after BC application: \n{}\n'.format(om_domain))
        n_om_domain = vorticity(om_domain,u_domain,v_domain,dx,dy,dt,nu,dim)
        logging.debug('Vorticity domain from transport equation: \n{}\n'.format(n_om_domain))
        si_domain = psi_calc(n_om_domain,si_domain,dx,beta,dim)
        logging.debug('StreamFunction domain from vorticity: \n{}\n'.format(si_domain))
        u_domain = u_calc(si_domain,u_domain,dy,dim)
        logging.debug('U velocity domain from streamfunction equation: \n{}\n'.format(u_domain))
        v_domain = v_calc(si_domain,v_domain,dx,dim)
        logging.debug('V velocity domain from streamfunction equation: \n{}\n'.format(v_domain))
        status = converge(om_domain,n_om_domain,criteria)
        np.copyto(om_domain,n_om_domain)
        iteration += 1

    logging.debug('Iterations completed. Status: {}\n'.format(status))

def discretize(n,m):
    return(np.zeros([n, m], dtype= float))

def vel_bound_condition(domain,bc,dim):
    m = dim[1]
    n = dim[0]
    for i in range(n):
        domain[i][0] = bc[0]
        domain[i][m-1] = bc[1]
        
    for j in range(m):
        domain[0][j] = bc[2]
        domain[n-1][j] = bc[3]
        
    return domain

#Central difference scheme
def om_calc(omega,u,v,dx,dy,dim,criteria):
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            omega[i][j] = ((v[i][j+1] - 2*v[i][j] + v[i][j-1])/(dx*dx)) - ((u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dy*dy))
    return omega



def vorticity_bound(omega,psi,u,v,dx,dy,dim):
    m = dim[1]
    n = dim[0]
    for j in range(1,n):
        omega[j][0] = 2*(psi[j][0]-psi[j][1])/(dx*dx) - 2*(v[j][0])/(dx)
        omega[j][m-1] = 2*(psi[j][m-1]-psi[j][m-2])/(dx*dx) + 2*(v[j][m-1])/(dx)

    for i in range(1,m):
        omega[0][i] = 2*(psi[0][i]-psi[1][i])/(dy*dy) + 2*(u[0][i])/(dy)
        omega[n-1][i] = 2*(psi[n-1][i]-psi[n-2][i])/(dy*dy) - 2*(u[n-1][i])/(dy)

    return omega

def vorticity(omega,u,v,dx,dy,dt,nu,dim):
    n_domain = discretize(dim[0],dim[1])
    m = dim[1]
    n = dim[0]
    for i in range(1,m-1):
        for j in range(1,n-1):
            n_domain[i][j] = omega[i][j] - (u[i][j]*(omega[i+1][j]-omega[i-1][j])/(2*dx))\
                - (v[i][j]*(omega[i][j+1]-omega[i][j-1])/(2*dy))\
                    + ((nu*dt)*(omega[i+1][j] - 2*omega[i][j] + omega[i-1][j])/(dx*dx))\
                        + ((nu*dt)*(omega[i][j+1] - 2*omega[i][j] + omega[i][j-1])/(dy*dy))
        
    return n_domain


#Elliptic equation (Poisson equation)
def psi_calc(omega,psi,dx,beta,dim):
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            psi[i][j] = ((psi[i][j+1]+psi[i][j-1]+\
                ((beta*beta)*(psi[i+1][j]+psi[i-1][j]))+\
                    (omega[i][j]*dx*dx))/\
                        (2*(1+beta*beta)))
    return psi

def u_calc(psi,u,dy,dim):    
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            u[i][j] = (psi[i+1][j]-psi[i-1][j])/(2*dy)
    return u

def v_calc(psi,v,dx,dim):    
    for i in range(1,dim[0]-1):
        for j in range(1,dim[1]-1):
            v[i][j] = (psi[i][j+1]-psi[i][j-1])/(2*dx)
    return v



def converge(o_omega,n_omega,criteria):
    convergence = np.subtract(o_omega,n_omega)
    if (np.amax(convergence,axis = None) <= criteria):
        return True
    else:
        return False

if (__name__ == "__main__"):
    prog_run()