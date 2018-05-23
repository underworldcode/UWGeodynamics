import numpy as np
import underworld as uw
import underworld.function as fn
import shapefile
import h5py
import numpy as np
import os
from .scaling import nonDimensionalize as nd
from .scaling import Dimensionalize
from .scaling import UnitRegistry as u
from ._swarm import Swarm
from scipy import spatial

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

        self.swarm = Swarm(mesh=mesh, particleEscape=particleEscape)
        self.swarm.add_particles_with_coordinates(points)
        self.advector = uw.systems.SwarmAdvector(swarm=self.swarm,
                                                 velocityField=velocityField,
                                                 order=2)
        indices = np.arange(self.swarm.particleLocalCount)
        rank = uw.rank()
        ranks = np.repeat(rank, self.swarm.particleLocalCount)
        pairs = np.array(list(zip(ranks, indices)), dtype=[("a", np.int32),
                                                           ("b", np.int32)])
        # Get rank
        self.global_index = self.swarm.add_variable(dataType="long", count=1)
        self.global_index.data[:, 0] = pairs.view(np.int64)

        self.tracked_field = list()

    def integrate(self, dt, **kwargs):
        """ Integrate swarm velocity in time """
        self.advector.integrate(dt, **kwargs)

    def add_tracked_field(self, value, name, units, dataType, count=1,
                          overwrite=True):
        """ Add a field to be tracked """
        if not isinstance(value, fn.Function):
            raise ValueError("%s is not an Underworld function")

        # Check that the tracer does not exist already
        for field in self.tracked_field:
            if (name == field["name"]) or (value == field["value"]):
                if not overwrite:
                    raise ValueError(""" %s name already exist or already tracked
                                     with a different name """ % name)
                else:
                    field["name"] = name
                    field["units"] = units
                    field["value"] = value
                    field["dataType"] = dataType
                    setattr(self, name, self.swarm.add_variable(dataType, count=count))
                    return

        self.tracked_field.append({"value":value,
                                   "name": name,
                                   "units": units,
                                   "dataType": dataType})
        setattr(self, name, self.swarm.add_variable(dataType, count=count))

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

    def save(self, outputDir, checkpointID, time):
        """ Save to h5 and create an xdmf file for each tracked field """

        # Save the swarm
        swarm_fname = self.name + '-%s.h5' % checkpointID
        swarm_fpath = os.path.join(outputDir, swarm_fname)
        sH = self.swarm.save(swarm_fpath, units=u.kilometers)

        filename = self.name + '-%s.xdmf' % checkpointID
        filename = os.path.join(outputDir, filename)

        # First write the XDMF header
        string = uw.utils._xdmfheader()
        string += uw.utils._swarmspacetimeschema(sH, swarm_fname, time.magnitude)

        # Save global index
        file_prefix = os.path.join(
            outputDir, self.name + '_global_index-%s' % checkpointID)
        handle = self.global_index.save('%s.h5' % file_prefix)
        string += uw.utils._swarmvarschema(handle, "global_index")

        # Save each tracked field
        for field in self.tracked_field:

            file_prefix = os.path.join(
                outputDir,
                self.name + "_" + field["name"] + '-%s' % checkpointID)

            obj = getattr(self, field["name"])
            obj.data[...] = field["value"].evaluate(self.swarm)
            handle = obj.save('%s.h5' % file_prefix, units=field["units"])

            # Add attribute to xdmf file
            string += uw.utils._swarmvarschema(handle, field["name"])

        # get swarm parameters - serially read from hdf5 file to get size
        h5f = h5py.File(name=swarm_fpath, mode="r")
        dset = h5f.get('data')
        if dset == None:
            raise RuntimeError("Can't find 'data' in file '{}'.\n".format(swarm_fname))
        globalCount = len(dset)
        dim = self.swarm.mesh.dim
        h5f.close()

        string += "\t<Attribute Type=\"Scalar\" Center=\"Node\" Name=\"Coordinates\">\n"
        string += """\t\t\t<DataItem Format=\"HDF\" NumberType=\" Float\"
                     Precision=\"8\" Dimensions=\"{0} {1}\">{2}:/data</DataItem>\n""".format(globalCount, dim, swarm_fname)
        string += "\t</Attribute>\n"

        # Write the footer to the xmf
        string += uw.utils._xdmffooter()

        # Write the string to file - only proc 0
        xdmfFH = open(filename, "w")
        xdmfFH.write(string)
        xdmfFH.close()


class Balanced_InflowOutflow(object):

    def __init__(self, vtop, top, pt1, pt2, ynodes=None,
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

        self.vtop = vtop
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


def circles_grid(radius, minCoord, maxCoord, npoints=72):


    if len(minCoord) == 2:
        # Create points on circle
        angles = np.linspace(0, 360, npoints)
        radius = nd(radius)
        x = radius * np.cos(np.radians(angles))
        y = radius * np.sin(np.radians(angles))

        # Calculate centroids
        xc = np.arange(nd(minCoord[0]), nd(maxCoord[0]) + radius, 2.*radius)
        yc = np.arange(nd(minCoord[1]) + radius, nd(maxCoord[1]), 2.*radius * np.sqrt(3)/2.)
        xc, yc = np.meshgrid(xc, yc)
        # Shift every other row by radius
        xc[::2,:] = xc[::2, :] + radius

        # Calculate coordinates of all circles points
        points = np.array(zip(xc.ravel(), yc.ravel()))[:, np.newaxis] + np.array(zip(x, y))
        x, y = points[:,:,0].ravel(), points[:,:,1].ravel()

        return x, y

    if len(minCoord) == 3:
        # Create points on circle
        theta = np.linspace(0, 180, npoints)
        phi = np.linspace(0, 360, npoints)
        radius = nd(radius)
        theta, phi = np.meshgrid(theta, phi)

        x = radius * np.sin(np.radians(theta.ravel())) * np.cos(np.radians(phi.ravel()))
        y = radius * np.sin(np.radians(theta.ravel())) * np.sin(np.radians(phi.ravel()))
        z = radius * np.cos(np.radians(theta.ravel()))

        # Calculate centroids
        xc = np.arange(nd(minCoord[0]) + radius, nd(maxCoord[0]) + radius, 2. * radius)
        yc = np.arange(nd(minCoord[1]) + radius, nd(maxCoord[1]) + radius, 2. * radius * np.sqrt(3)/2.)
        zc = np.arange(nd(minCoord[2]) + radius, nd(maxCoord[2]) + radius, 2. * radius * np.sqrt(3)/2.)
        xc, yc, zc = np.meshgrid(xc, yc, zc)
        # Shift every other row by radius
        yc[:,::2,:] += radius
        zc[::2,:,:] += radius

        # Calculate coordinates of all circles points
        points = np.array(zip(xc.ravel(), yc.ravel(), zc.ravel()))[:, np.newaxis] + np.array(zip(x, y, z))
        x = points[:,:,0].ravel()
        y = points[:,:,1].ravel()
        z = points[:,:,2].ravel()


        return x, y, z

def circle_points_tracers(radius, centre=tuple([0., 0.]), npoints=72):
    angles = np.linspace(0, 360, npoints)
    radius = nd(radius)
    x = radius * np.cos(np.radians(angles)) + nd(centre[0])
    y = radius * np.sin(np.radians(angles)) + nd(centre[1])
    return x, y

def sphere_points_tracers(radius, centre=tuple([0., 0., 0.]), npoints=30):
    theta = np.linspace(0, 180, npoints)
    phi = np.linspace(0, 360, npoints)
    radius = nd(radius)
    theta, phi = np.meshgrid(theta, phi)

    x = radius * np.sin(np.radians(theta.ravel())) * np.cos(np.radians(phi.ravel()))
    y = radius * np.sin(np.radians(theta.ravel())) * np.sin(np.radians(phi.ravel()))
    z = radius * np.cos(np.radians(theta.ravel()))

    x += nd(centre[0])
    y += nd(centre[1])
    z += nd(centre[2])

    return x, y, z


class Nearest_neigbhors_projector(object):

    def __init__(self, mesh, swarm, swarm_variable, mesh_variable, dtype):
        self.mesh = mesh
        self.swarm = swarm
        self.swarm_variable = swarm_variable
        self.mesh_variable = mesh_variable

    def solve(self):
        tree = spatial.KDTree(self.swarm.particleCoordinates.data)
        ids = tree.query(self.mesh.data)
        pts = self.swarm.particleCoordinates.data[ids, :]
        self.mesh_variable.data[...] = self.swarm_variable.evaluate(pts)


def fn_Tukey_window(r, centre, width, top, bottom):
    """ Define a tuckey window

    A Tukey window is a rectangular window with the first and last r/2
    percent of the width equal to parts of a cosine.
    see tappered cosine function

    """

    centre = nd(centre)
    width = nd(width)
    top = nd(top)
    bottom = nd(bottom)

    x = fn.input()[0]
    y = fn.input()[1]

    start = centre - 0.5 * width
    xx = (x - start) / width
    x_conditions = [((0. <= xx) & (xx < r /2.0), 0.5 * (1.0 + fn.math.cos(2. * np.pi / r *(xx - r / 2.0)))),
                    ((r / 2.0 <= xx) & (xx < 1.0 - r / 2.0), 1.0),
                    ((1.0 - r / 2.0 <= xx) & (xx <= 1.0), 0.5 * (1. + fn.math.cos(2. * np.pi / r *(xx + r / 2.0)))),
                    (True, 0.0)]

    x_conditions =  fn.branching.conditional(x_conditions)

    y_conditions = fn.branching.conditional([((y >= bottom) & (y <= top), 1.0), (True, 0.0)])
    return x_conditions * y_conditions


import pylab as plt

class NonLinearBlock(object):
    def __init__(self, string):
        self.string = string
        self.data = dict()
        self.data["Pressure Solve times"] = self.get_vals(["Pressure Solve"], 3)
        self.data["Final V Solve times"] = self.get_vals(["Final V Solve"], 4)
        self.data["Total BSSCR times"] = self.get_vals(["Total BSSCR Linear solve time"], 5)
        self.data["Residuals"] = self.get_vals(["Iteration", "Residual", "Tolerance"], 9)
        self.data["Iterations"] = self.get_vals(["Non linear solver - iteration"], -1, func=int)
        self.data["Solution Time"] = self.get_vals(["solution time"], 5)


    def get_vals(self, FINDSTRING, pos, func=float):
        f = self.string.splitlines()
        vals = [func(line.split()[pos]) for line in f if all([F in line for F in FINDSTRING])]
        return vals

class LogFile(object):
    def __init__(self, filename):
        self.filename = filename
        self.nonLinear_blocks = self.get_nonLinear_blocks()
        self.pressure_solve_times = list()
        for obj in self.nonLinear_blocks:
            self.pressure_solve_times += obj.data["Pressure Solve times"]
        self.finalV_solve_times = list()
        for obj in self.nonLinear_blocks:
            self.finalV_solve_times += obj.data["Final V Solve times"]
        self.total_BSSCR_times = list()
        for obj in self.nonLinear_blocks:
            self.total_BSSCR_times += obj.data["Total BSSCR times"]
        self.residuals = list()
        for obj in self.nonLinear_blocks:
            self.residuals += obj.data["Residuals"]
        self.iterations = list()
        for obj in self.nonLinear_blocks:
            self.iterations.append(obj.data["Iterations"][-1])
        self.solution_times = list()
        for obj in self.nonLinear_blocks:
            self.solution_times += obj.data["Solution Time"]

    def get_nonLinear_blocks(self):
        non_linear_blocks = list()
        with open(self.filename, "r") as f:
            block = ""
            inBlock = False
            step=0
            for line in f:
                if "Non linear solver" in line:
                    inBlock = True
                    step += 1
                if inBlock:
                    if "Converged" not in line:
                        block += line
                    else:
                        block += line
                        inBlock = False
                        block = NonLinearBlock(block)
                        non_linear_blocks.append(block)
                        block=""
        self.nonLinear_blocks = non_linear_blocks
        return self.nonLinear_blocks
