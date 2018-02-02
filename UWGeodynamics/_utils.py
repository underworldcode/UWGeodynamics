import numpy as np
import underworld as uw
from .scaling import nonDimensionalize as nd

class PressureSmoother(object):

    def __init__(self, mesh, pressureField):

        self.mesh = mesh
        self.pressureField = pressureField

        self.NodePressure = uw.mesh.MeshVariable(self.mesh, nodeDofCount=1)

        self.Cell2Nodes = uw.utils.MeshVariable_Projection(
            self.NodePressure, self.pressureField, type=0)
        self.Nodes2Cell = uw.utils.MeshVariable_Projection(
            self.pressureField, self.NodePressure, type=0)

    def smooth(self):

        self.Cell2Nodes.solve()
        self.Nodes2Cell.solve()


class PassiveTracers(object):

    def __init__(self, mesh, velocityField, name=None, vertices=None,
                 particleEscape=True):

        self.name = name
        self.particleEscape = particleEscape

        for dim in range(len(vertices)):
            vertices[dim] = nd(vertices[dim])

        sizes = np.array([np.array(x).size for x in vertices])
        points = np.zeros((sizes.max(), len(vertices)))

        for dim in range(len(vertices)):
            points[:, dim] = vertices[dim]

        self.swarm = uw.swarm.Swarm(mesh=mesh, particleEscape=particleEscape)
        self.swarm.add_particles_with_coordinates(points)
        self.advector = uw.systems.SwarmAdvector(swarm=self.swarm,
                                                 velocityField=velocityField,
                                                 order=2)

    def integrate(self, dt, **kwargs):

        self.advector.integrate(dt, **kwargs)


class Balanced_InflowOutflow(object):

    def __init__(self, Vtop, top, pt1, pt2, ynodes=None, 
                 tol=1e-12, nitmax=200, 
                 nitmin=3, default_vel=0.0):
        """ Calculate Bottom velocity as for Huismans et al. velocity boundary
            conditions such as the total volume is conserved.
            Return the velocity profile.

            NOTE. The current version implies uniform dy.

            Input:

            Vtop     Top Velocity condition
            pt1, pt2 Top and Bottom location of the transition zone where the
                     velocity linearly decreases from Vtop to VBottom.
            ynodes   Coordinates of the nodes in the y-direction.

            Output:

            velocities numpy array.

            """

        self.vtop = Vtop
        self.top = top
        self.pt1 = pt1
        self.pt2 = pt2
        self.ynodes = ynodes
        self.tol = tol
        self.nitmax = nitmax
        self.nitmin = nitmin
        self.default_vel = default_vel


    def _get_side_flow(self):
   
        Vtop = nd(self.vtop)
        top = nd(self.top)
        pt1 = nd(self.pt1)
        pt2 = nd(self.pt2)
        y = nd(self.ynodes)
        tol = self.tol
        nitmax = self.nitmax
        nitmin = self.nitmin
        default_vel = nd(self.default_vel)

        # locate index of closest node coordinates
        top_idx = np.argmin((y - top)**2)
        pt1_idx = np.argmin((y - pt1)**2)
        pt2_idx = np.argmin((y - pt2)**2)
    
        # do some initialization
        velocity = np.ones(y.shape) * Vtop
        Vmin = -Vtop
        Vmax = 0.0
        N = 0
    
        dy = np.diff(y)
        budget = 0.5 * (velocity[1:] + velocity[:-1]) * dy
        prev = np.copy(budget)
    
        # The following loop uses a kind of bissection approach
        # to look for the suitable value of Vbot.
        while True:
    
            Vbot = (Vmin + Vmax) / 2.0
    
            for i in range(len(y)):
                if i > top_idx:
                    velocity[i] = 0.0
                if i >= pt1_idx and i <= top_idx:
                    velocity[i] = Vtop
                if i <= pt2_idx:
                    velocity[i] = Vbot
                if i < pt1_idx and i > pt2_idx:
                    velocity[i] = (Vtop - Vbot) / (y[pt1_idx] - y[pt2_idx]) * (y[i] - y[pt2_idx]) + Vbot
    
            budget = 0.5 * (velocity[1:] + velocity[:-1]) * dy
    
            if np.abs(np.sum(budget) - np.sum(prev)) < tol and N > nitmin:
                velocity[top_idx + 1:] = default_vel
                self.budget = np.sum(budget)
                return velocity
            else:
                ctol = np.abs(np.sum(budget) - np.sum(prev))
                N += 1
                prev = np.copy(budget)
    
            if Vtop < 0.0:
                if np.sum(budget) < 0.0:
                    Vmax = Vbot
                else:
                    Vmin = Vbot
            else:
                if np.sum(budget) > 0.0:
                    Vmax = Vbot
                else:
                    Vmin = Vbot
    
        velocity[top_idx + 1:] = default_vel

        self.budget = np.sum(budget)
        
        return velocity
    
