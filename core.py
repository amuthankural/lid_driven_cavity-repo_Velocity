import sys
import logging
import json
import numpy as np


import velocity
import vorticity
import streamfunc
import converge


np.set_printoptions(precision=2)
np.set_printoptions(threshold=sys.maxsize)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')
file_handler = logging.FileHandler('ftcs.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

def prog_run():
    with open('./input.json','r') as f:
        input           = json.load(f)

    """-----------Imported data----------"""
    l_i                 = input["lx_i"]                 #Length of the domain
    b_j                 = input["by_j"]                 #Breadth of the domain
    m_i                 = input["m_i"]                  #Number of grid points along 'x'
    n_j                 = input["n_j"]                  #Number of grid points along 'y'
    Re                  = input["Re"]                   #Reynolds Number
    bc_u                = [input["u_L"],input["u_R"],input["u_B"],input["u_T"]]
    bc_v                = [input["v_L"],input["v_R"],input["v_B"],input["v_T"]]
    criteria            = input["criteria"]
    nu                  = input["nu"]
    dt                  = input["dt"]
    status              = False
    iteration           = 1
    logger.info('Data imported from JSON file\n---------------------------------------------------------------------')
    logger.info('\nDomain data:\n\tDomain length:\t{}\n\tDomain breadth:\t{}\n\tX divisions:\t{}\n\tY divisions:\t{}'.format(l_i,b_j,m_i,n_j))
    logger.info('\n---------------------------------------------------------------------')
    logger.info('Flow properties:')
    logger.info('\nu Boundary conditions:\n\tLeft wall:\t{}\n\tRight wall:\t{}\n\tBottom wall:\t{}\n\tTop wall:\t{}'\
        .format(bc_u[0],bc_u[1],bc_u[2],bc_u[3]))
    logger.info('\nv Boundary conditions:\n\tLeft wall:\t{}\n\tRight wall:\t{}\n\tBottom wall:\t{}\n\tTop wall:\t{}'\
        .format(bc_v[0],bc_v[1],bc_v[2],bc_v[3]))
    logger.info('Reynold\'s number:\t{}'.format(Re))
    logger.info('Criteria for convergence:\t{}\n---------------------------------------------------------------------'.format(criteria))
    logger.info('Fluid properties:')
    logger.info('Viscosity:\t{}'.format(nu))


    """-----------Calculated data----------"""
    dx                  = l_i/(m_i-1)                   #Grid length along 'x' direction
    dy                  = b_j/(n_j-1)                   #Grid length along 'y' direction
    beta                = dx/dy
    
    logger.info('dx Value:\t{}'.format(dx))
    logger.info('dy Value:\t{}'.format(dy))


    """---------------Domains--------------"""
    n_om_domain             = discretize(n_j,m_i)
    o_om_domain             = discretize(n_j,m_i)
    n_psi_domain            = discretize(n_j,m_i)
    o_psi_domain            = discretize(n_j,m_i)
    u_domain                = discretize(n_j,m_i)
    v_domain                = discretize(n_j,m_i)
    logger.info('Domain arrays created\n---------------------------------------------------------------------')


    u_domain =velocity.boundary(u_domain,bc_u,m_i,n_j)
    logger.info('U Velocity domain:\n{}'.format(u_domain))
    v_domain = velocity.boundary(v_domain,bc_v,m_i,n_j)
    logger.info('V Velocity domain:\n{}'.format(v_domain))
    while iteration <= 5:
        n_om_domain = vorticity.transposrt(n_om_domain,o_om_domain,o_psi_domain,dx,dy,dt,m_i,n_j,nu)
        logger.info('Vorticity domain from transport equation:\n{}'.format(n_om_domain))
        n_om_domain = vorticity.boundary(n_om_domain,o_psi_domain,u_domain,v_domain,dx,dy,m_i,n_j)
        logger.info('Vorticity domain Boundary condition:\n{}'.format(n_om_domain))
        n_psi_domain = streamfunc.poisson(n_om_domain,o_psi_domain,dx,m_i,n_j)
        logger.info('Stream function from poisson equation:\n{}'.format(n_psi_domain))
        status = converge.check(o_om_domain,n_om_domain,criteria)
        logger.info('Convergence status for iteration number {} is:\n{}'.format(iteration,status))
        iteration += 1
        logger.info('\n==============================================================================')




def discretize(n,m):
    return(np.zeros([n, m], dtype= float))


if (__name__ == "__main__"):
    prog_run()