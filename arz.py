#!/usr/bin/env python
# encoding: utf-8
import numpy as np

tau = 5.

def hesitation(rho):
    return 25. * rho**0.2/(1.-rho)**0.1

def g(y):
    return np.sqrt(1. + (10.*(y-1./3))**2)

def desired_velocity(rho):
    return 1.4976*(g(0.) + (g(1.)-g(0.))*rho - g(rho))/rho

def step_velocity_relaxation(solver,state,dt):
    """Compute velocity relaxation term in conserved variables:

        $$(\rho (U(\rho) + h(\rho))-q)/tau$$
    """
    rho = state.q[0,:]
    q   = state.q[1,:]

    velocity = q/rho - hesitation(rho)
    src = np.zeros(state.q.shape)
    state.q[1,:] = state.q[1,:] + dt*(desired_velocity(rho) - velocity)/tau
    return src

def dq_velocity_relaxation(solver,state,dt):
    """Compute velocity relaxation term in conserved variables:

        $$(\rho (U(\rho) + h(\rho))-q)/tau$$
    """
    rho = state.q[0,:]
    q   = state.q[1,:]

    velocity = q/rho - hesitation(rho)
    src = np.zeros(state.q.shape)
    src[1,:] = dt*(desired_velocity(rho) - velocity)/tau
    return src


def arz(use_petsc=False,solver_type='sharpclaw',iplot=False,htmlplot=False,outdir='./_output',weno_order=5,rp_type='fwave',IC='wiggles',num_cells=1000):
    """
    This example solves the 1-dimensional AW-Rascle traffic model.
    The conserved variables are:
    q[0,:]: \rho
    q[0,:]: \rho(u + h(\rho))

    where $\rho$ is the density and $h(\rho)$ is the hesitation function.
    """

    #=================================================================
    # Import the appropriate classes, depending on the options passed
    #=================================================================
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D()
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.step_source=step_velocity_relaxation
        solver.order=2
        solver.cfl_max = 0.45
        solver.cfl_desired = 0.4
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
        solver.weno_order=weno_order
        solver.time_integrator = 'SSP104'
        solver.cfl_max = 2.45
        solver.cfl_desired = 2.4
        solver.dq_src=dq_velocity_relaxation
    else: raise Exception('Unrecognized value of solver_type.')

    solver.num_waves=2
    solver.num_eqn=2

    if rp_type == 'fwave':
        import rp1_arz_traffic
        solver.rp = rp1_arz_traffic
        solver.fwave = True
    elif rp_type == 'hll':
        import rp1_hll
        solver.rp = rp1_hll
        solver.fwave = False

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    x = pyclaw.Dimension(0.0,500.0,num_cells)
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)
    claw = pyclaw.Controller()

    #========================================================================
    # Set the initial condition
    #========================================================================
    xc=domain.grid.x.centers
    if IC=='wiggles':
        state.q[0,:] = 0.6 + 0.005*np.sin(2*np.pi*xc/500.) + 0.005*np.sin(24*np.pi*xc/500.)
        state.q[1,:] = state.q[0,:]*hesitation(state.q[0,:]) # Zero velocity
        claw.tfinal = 600.
        claw.num_output_times = 180
    elif IC=='rp1':
        state.q[0,:] = 0.9*(xc<250.)+0.1*(xc>250.)
        state.q[1,:] = state.q[0,:]*(1.+hesitation(state.q[0,:]))
        solver.dq_src=None
        solver.bc_lower[0] = pyclaw.BC.extrap
        solver.bc_upper[0] = pyclaw.BC.extrap
        claw.tfinal = 1.
        claw.num_output_times = 10


    #========================================================================
    # Set up the controller object
    #========================================================================
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir

    # Plot results
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(arz)
