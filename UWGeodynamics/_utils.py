import numpy as np
import underworld as uw
import shapefile
import numpy as np
import os
from .scaling import nonDimensionalize as nd
from .scaling import Dimensionalize
from .scaling import UnitRegistry as u

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
                 particleEscape=True, shapeType="line"):

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


    def write_to_shapefile(self, filename, units=None, overwrite=False):

        if os.path.exists(filename) and not overwrite:
            r = shapefile.Reader(filename)
            w = shapefile.Writer(r.shapeType)
            # Copy over the existing dbf fields
            w.fields = list(r.fields)
            # Copy over the existing dbf records
            w.records.extend(r.records())
            # Copy over the existing polygons
            w._shapes.extend(r.shapes())
 
        else:

            w = shapefile.Writer(shapeType=shapefile.POLYLINEZ)
            w.field("name", "C")
            w.field("units", "C")

        fact = 1.0
        if units:
            fact = Dimensionalize(1.0, units=units)
            fact = fact.magnitude

        x = self.swarm.particleCoordinates.data[:,0] * fact
        y = self.swarm.particleCoordinates.data[:,1] * fact
        line = zip(x,y)
        w.poly(parts=[line])
        w.record(self.name, str(units))
        w.save(filename)


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

    def to_json(self):
        d = {}
        d["type"] = "Balanced_InflowOutflow"
        d["vtop"] = str(self.vtop)
        d["top"] = str(self.top)
        d["pt1"] = str(self.pt1)
        d["pt2"] = str(self.pt2)
        d["ynodes"] = self.ynodes
        d["tol"] = str(self.tol)
        d["nitmax"] = str(self.nitmax)
        d["mitmin"] = str(self.nitmin)
        d["default_vel"] = str(self.default_vel)
        return d


class MoveImporter(object):

    def __init__(self, filename, units):

        self.filename = filename
        self.shapes = self.records()
        self.units = units
        self.attributes = None
        self.extent = None
        self.names = []

        xmin = float('inf') * units
        xmax = -float('inf') * units
        ymin = float('inf') * units
        ymax = -float('inf') * units

        _, self.filetype = os.path.splitext(filename)

        for record in self.records():
            self.names.append(record["properties"]["Name"])
            if record["extent"][0][0].magnitude < xmin.magnitude:
                xmin = record["extent"][0][0]
            if record["extent"][0][1].magnitude > xmax.magnitude:
                xmax = record["extent"][0][1]
            if record["extent"][1][0].magnitude < ymin.magnitude:
                ymin = record["extent"][1][0]
            if record["extent"][1][1].magnitude > ymax.magnitude:
                ymax = record["extent"][1][1]

        self.names = np.unique(self.names)

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        self.coords = [(xmin, ymin), (xmax, ymax)]

        self.generator = self.records()

    def records(self):

        if self.units:
            units = self.units
        else:
            units = u.dimensionless

        reader = shapefile.Reader(self.filename)
        fields = reader.fields[1:]
        field_names = [field[0] for field in fields]

        for sr in reader.shapeRecords():
            atr = dict(zip(field_names, sr.record))
            coords = np.array(sr.shape.points)
            coords[:,-1] = sr.shape.z
            xextent = (coords[:,0].min(), coords[:,0].max())
            yextent = (coords[:,1].min(), coords[:,1].max())
            yield dict(coordinates=coords * units, properties=atr,
                       extent=[xextent * units, yextent * units])



