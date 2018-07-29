import abc
try:
    from .linkage import SPM
except:
    pass
import underworld.function as fn
from scaling import nonDimensionalize as nd

import numpy as np
from scipy.linalg import solve as linSolve
import scipy.signal as sig
from scipy.interpolate import interp1d
from mpi4py import MPI
comm = MPI.COMM_WORLD
CPUsize = comm.Get_size()

ABC = abc.ABCMeta('ABC', (object,), {})

class SurfaceProcesses(ABC):

    def __init__(self, Model=None):

        self.Model = Model

    @property
    def Model(self):
        return self._Model

    @Model.setter
    def Model(self, value):
        self._Model = value
        if value:
            self._init_model()

    @abc.abstractmethod
    def _init_model(self):
        pass

    @abc.abstractmethod
    def solve(self, dt):
        pass


class Badlands(SurfaceProcesses):
    """ A wrapper class for Badlands Linkage"""

    def __init__(self, airIndex,
                 sedimentIndex, XML, resolution, checkpoint_interval,
                 surfElevation=0., verbose=True, Model=None, restartFolder=None,
                 restartStep=None, timeField=None, surfaceTracers=None):

        self.airIndex = airIndex
        self.sedimentIndex = sedimentIndex
        self.XML = XML
        self.resolution = resolution
        self.checkpoint_interval = checkpoint_interval
        self.surfElevation = surfElevation
        self.verbose = verbose
        self.restartFolder = restartFolder
        self.restartStep = restartStep
        self.surfaceTracers = surfaceTracers

        self.timeField = timeField
        self.Model = Model

    def _init_model(self):
        self.mesh = self._Model.mesh
        self.velocityField = self._Model.velocityField
        self.swarm = self._Model.swarm
        self.materialField = self._Model.materialField

        self._BadlandsModel = SPM(self.mesh, self.velocityField, self.swarm,
                                 self.materialField, self.airIndex, self.sedimentIndex,
                                 self.XML, nd(self.resolution),
                                 nd(self.checkpoint_interval),
                                 nd(self.surfElevation), self.verbose,
                                 self.restartFolder, self.restartStep)
        return

    def solve(self, dt):
        self._BadlandsModel.solve(dt)
        return


class ErosionThreshold(SurfaceProcesses):

    def __init__(self, air=None, threshold=None, surfaceTracers=None, Model=None):

        super(ErosionThreshold, self).__init__(Model=Model)

        self.Model = Model
        self.threshold = threshold
        self.air = air
        self.surfaceTracers=surfaceTracers
        self.Model = Model

    def _init_model(self):

        materialField = self.Model.materialField

        materialMap = {}
        for material in self.air:
            materialMap[material.index] = 1.0

        isAirMaterial = fn.branching.map(fn_key=materialField, mapping=materialMap, fn_default=0.0)

        belowthreshold = [(((isAirMaterial < 0.5) & (fn.input()[1] > nd(self.threshold))), self.air[0].index),
                         (True, materialField)]

        self._fn = fn.branching.conditional(belowthreshold)

    def solve(self, dt):

        if not self.Model:
            raise ValueError("Model is not defined")

        self.Model.materialField.data[:] = self._fn.evaluate(self.Model.swarm)
        if self.surfaceTracers:
            if self.surfaceTracers.swarm.particleCoordinates.data.size > 0:
                coords = self.surfaceTracers.swarm.particleCoordinates
                coords.data[coords.data[:,-1] > nd(self.threshold), -1] = nd(self.threshold)
        return


class SedimentationThreshold(SurfaceProcesses):

    def __init__(self, air=None, sediment=None,
                 threshold=None, timeField=None, Model=None,
                 surfaceTracers=None):

        super(SedimentationThreshold, self).__init__(Model=Model)

        self.timeField = timeField
        self.air = air
        self.sediment = sediment
        self.threshold = threshold
        self.surfaceTracers = surfaceTracers
        self.Model = Model

    def _init_model(self):

        materialField = self.Model.materialField

        materialMap = {}
        for material in self.air:
            materialMap[material.index] = 1.0

        isAirMaterial = fn.branching.map(fn_key=materialField, mapping=materialMap, fn_default=0.0)

        belowthreshold = [(((isAirMaterial > 0.5) & (fn.input()[1] < nd(self.threshold))), 0.),
                         (True, 1.)]

        self._change_material = fn.branching.conditional(belowthreshold)

        conditions = [(self._change_material < 0.5, self.sediment[0].index),
                      (True, materialField)]

        self._fn = fn.branching.conditional(conditions)

    def solve(self, dt):

        if not self.Model:
            raise ValueError("Model is not defined")

        if self.timeField:
            fn = self._change_material * self.timeField
            self.timeField.data[...] = fn.evaluate(self.Model.swarm)

        self.Model.materialField.data[:] = self._fn.evaluate(self.Model.swarm)

        if self.surfaceTracers:
            if self.surfaceTracers.swarm.particleCoordinates.data.size > 0:
                coords = self.surfaceTracers.swarm.particleCoordinates
                coords.data[coords.data[:,-1] < nd(self.threshold), -1] = nd(self.threshold)


class ErosionAndSedimentationThreshold(SedimentationThreshold, ErosionThreshold):

    def __init__(self, air=None, sediment=None,
                 threshold=None, timeField=None,
                 surfaceTracers=None, Model=None):

        super(ErosionAndSedimentationThreshold, self).__init__(Model=Model)

        self.timeField = timeField
        self.air = air
        self.sediment = sediment
        self.threshold = threshold
        self.surfaceTracers = surfaceTracers
        self.Model = Model

    def _init_model(self):

        ErosionThreshold._init_model(self)
        SedimentationThreshold._init_model(self)

    def solve(self, dt):

        ErosionThreshold.solve(self, dt)
        SedimentationThreshold.solve(self, dt)

class BasicHillSlopeDiffsuion2d(object):
    def __init__(self,
                 Model=None,
                 airIndex=None,
                 sedimentIndex=None,
                 diffusivity=None,
                 interfaceHeight=0.,
                 timeField=None,
                 filterTopo=False,
                 verbose=True):

        # Create references to Model variables
        self.materialField = Model.materialField
        self.airIndex = airIndex
        self.sedimentIndex = sedimentIndex
        self.Ks = diffusivity
        self.mesh = Model.mesh
        self.velocityField = Model.velocityField
        self.swarm = Model.swarm
        self.minX = nd(Model.minCoord[0])
        self.maxX = nd(Model.maxCoord[0])

        # Define the number of topographic markers, 4 times the Model mesh resolution
        self.topoNum = 4 * Model.elementRes[0] + 1
        self.topostp = (self.maxX - self.minX) / (self.topoNum - 1)

        # initiate the 1d FCT Grid for topographic diffusion
        self.gridt = np.zeros((6, self.topoNum))
        self.gridt[0, :] = np.linspace(self.minX, self.maxX, self.topoNum)
        self.gridt[1, :] = nd(interfaceHeight)
        self.filterTopo = filterTopo
        self.verbose = verbose

    def SurfaceVeloEval(self, mesh=None, velocityField=None):

        minX = self.minX
        maxX = self.maxX

        self.gridt[3:6, :] = 0.0

        tmp = np.where(
            (self.gridt[0, :] >= minX) & (self.gridt[0, :] <= maxX) &
            (self.gridt[0, :] >= mesh.data[0:mesh.nodesLocal, 0].min()) &
            (self.gridt[0, :] <= mesh.data[0:mesh.nodesLocal, 0].max())
            & (self.gridt[1, :] <= mesh.data[0:mesh.nodesLocal, 1].max()))[0]

        if len(tmp) > 0:
            tmp2 = velocityField.evaluate(np.squeeze(self.gridt[0:2, tmp]).T)
            self.gridt[3, tmp] = tmp2.T[0, :]
            self.gridt[4, tmp] = tmp2.T[1, :]

            tmp = np.where(
                (self.gridt[0, :] > minX) & (self.gridt[0, :] < maxX) &
                ((self.gridt[0, :] == mesh.data[0:mesh.nodesLocal, 0].min())
                 | (self.gridt[0, :] == mesh.data[0:mesh.nodesLocal, 0].max()))
            )[0]
            # boundary between two cpus, there velocity is reduced
            if len(tmp) > 0:
                # import ipdb; ipdb.set_trace()
                print 'hgn', tmp, self.gridt[0:2, tmp], np.squeeze(
                    self.gridt[0:2, tmp]).T
                if len(tmp) == 1:
                    tmp2 = velocityField.evaluate((self.gridt[0, tmp][0],
                                                   self.gridt[1, tmp][0]))
                else:
                    tmp2 = velocityField.evaluate(
                        np.squeeze(self.gridt[0:2, tmp]).T)
                self.gridt[3, tmp] = tmp2.T[0, :] / 2.
                self.gridt[4, tmp] = tmp2.T[1, :] / 2.

    def SurfaceProcess(self, dt):

        Ks = self.Ks
        topoNum = self.topoNum
        topostp = self.topostp
        minX = self.minX
        maxX = self.maxX
        # refer to Collision.m in Chapter_17 of Gerya_numerical_geodynamics book
        # first advect topography vertically
        # and diffuse topography (downhill diffusion)
        L = np.zeros((topoNum, topoNum))
        R = np.zeros((topoNum, 1))
        # first point: symmetry
        L[0, 0] = 1.
        L[0, 1] = -1.
        R[0] = 0.0
        # from IPython.core.debugger import Tracer; Tracer()()
        # Intermediate Points
        K2 = Ks * dt / topostp**2
        for i1 in range(1, topoNum - 1):
            # Internal points
            if (self.gridt[0, i1] >= minX and self.gridt[0, i1] <= maxX):
                L[i1, i1 - 1] = -K2
                L[i1, i1] = 1 + 2 * K2
                L[i1, i1 + 1] = -K2
                R[i1] = self.gridt[1, i1] + self.gridt[4, i1] * dt
            else:
                # left of the left boundary
                if (self.gridt[0, i1] < minX):
                    L[i1, i1] = 1.
                    L[i1, i1 + 1] = -1.
                    R[i1] = 0

                # right of the right boundary
                if (self.gridt[0, i1] > maxX):
                    L[i1, i1] = 1.
                    L[i1, i1 - 1] = -1.
                    R[i1] = 0

        # last point: symmetry
        L[topoNum - 1, topoNum - 1] = 1.
        L[topoNum - 1, topoNum - 2] = -1.
        R[topoNum - 1] = 0.

        # solve matrix
        self.gridt[1, :] = np.squeeze(linSolve(L, R))
        # Second, advect topography horizontally
        vxmax = max(np.abs(self.gridt[
            3, :]))  # maximum horizontal velocity at topography nodes
        # defining topography advection timestep
        ntSP = 1
        dtSP = dt
        if vxmax > 0:
            dtSP = min(topostp / vxmax, dt)
            if dtSP < dt:
                ntSP = np.ceil(dt / dtSP)
                dtSP = dt / ntSP

        # define FCT parameter MU
        mu = 1.0 / 8
        # advect topography with FCT
        for i1 in range(ntSP):
            # step 0: set new profile
            self.gridt[2, :] = self.gridt[1, :]
            # step 1: Transport + numerical diffusion stage
            for i2 in range(1, topoNum - 1):
                # define FCT parameters EPS and NU
                eps = self.gridt[3, i2] * dtSP / topostp
                nu = 1. / 8 + eps**2 / 2.
                # change topo
                self.gridt[2, i2] = self.gridt[1, i2] - eps / 2 * (
                    self.gridt[1, i2 + 1] - self.gridt[1, i2 - 1]) + nu * (
                        self.gridt[1, i2 + 1] - 2 * self.gridt[1, i2] +
                        self.gridt[1, i2 - 1])

            # step 2: anti-difussion stage
            # anti-diffusion flow for the first cell
            self.gridt[5, 0] = 0
            for i2 in range(1, topoNum - 2):
                # corrected antidiffusion flow for current cell
                delt0 = self.gridt[2, i2] - self.gridt[2, i2 - 1]
                delt1 = self.gridt[2, i2 + 1] - self.gridt[2, i2]
                delt2 = self.gridt[2, i2 + 2] - self.gridt[2, i2 + 1]
                s = np.copysign(1.0, delt1)
                self.gridt[5, i2] = s * max(
                    0.0, min(min(s * delt2, s * delt0), mu * abs(delt1)))
                self.gridt[
                    1,
                    i2] = self.gridt[2, i2] - self.gridt[5,
                                                         i2] + self.gridt[5, i2
                                                                          - 1]

        # Filter/Moving average to remove smale scale instabilities
        # for certain values of Ks or when dt is large

    #
        if self.filterTopo:
            self.gridt[1, :] = sig.savgol_filter(
                self.gridt[1, :], 3, 1, mode='nearest')
        return

    def ErosionAndSedimentation(self):

        airIndex = self.airIndex
        sedimentIndex = self.sedimentIndex

        # generate an interpolation function, nearest seems to be the fastest option, refer to linkage module.
        surface_function = interp1d(
            self.gridt[0, :], self.gridt[1, :], kind='nearest')
        swarm_coords = self.swarm.particleCoordinates.data
        surface_ycoord = surface_function(swarm_coords[:, 0])
        material_flags = swarm_coords[:, 1] < surface_ycoord

        # convert air to sediment
        sedimented_mask = np.logical_and(
            np.in1d(self.materialField.data, airIndex), material_flags)
        self.materialField.data[sedimented_mask] = sedimentIndex

        # convert sediment to air
        eroded_mask = np.logical_and(
            ~np.in1d(self.materialField.data, airIndex), ~material_flags)
        self.materialField.data[eroded_mask] = airIndex

        return

    def solve(self, dt):
        if comm.rank == 0 and self.verbose:
            purple = "\033[0;35m"
            endcol = "\033[00m"
            print(purple +
                  "Processing surface with BasicHillSlopeDiffsuion2d" + endcol)

        self.SurfaceVeloEval(mesh=self.mesh, velocityField=self.velocityField)
        self.gridt[3:5, :] = comm.allreduce(self.gridt[3:5, :], op=MPI.SUM)

        comm.barrier()
        if comm.rank == 0:
            self.SurfaceProcess(dt)
        self.gridt[1, :] = comm.bcast(self.gridt[1, :], root=0)
        comm.barrier()
        self.ErosionAndSedimentation()
        comm.barrier()

        if comm.rank == 0 and self.verbose:
            purple = "\033[0;35m"
            endcol = "\033[00m"
            print(purple +
                  "Processing surface with BasicHillSlopeDiffsuion2d" + endcol)
        return
