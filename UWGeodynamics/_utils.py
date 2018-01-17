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

        points = np.zeros((len(vertices[0]), len(vertices)))

        for dim in range(len(vertices)):
            points[:, dim] = vertices[dim]

        self.swarm = uw.swarm.Swarm(mesh=mesh, particleEscape=particleEscape)
        self.swarm.add_particles_with_coordinates(points)
        self.advector = uw.systems.SwarmAdvector(swarm=self.swarm,
                                                 velocityField=velocityField,
                                                 order=2)

    def integrate(self, dt, **kwargs):

        self.advector.integrate(dt, **kwargs)
