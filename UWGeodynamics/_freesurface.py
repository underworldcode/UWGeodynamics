import underworld as uw
from UWGeodynamics import nd
from scipy.interpolate import interp1d

class FreeSurfaceProcessor(object):
    """FreeSurfaceProcessor"""

    def __init__(self, Model):
        """__init__

        Parameters
        ----------

        Model : UWGeodynamics Model

        """
        self.Model = Model

        # Create the tools
        self.TField = self.Model.mesh.add_variable( nodeDofCount=1 )
        self.TField.data[:, 0] = self.Model.mesh.data[:, 1]

        self.top = self.Model.top_wall
        self.bottom = self.Model.bottom_wall

        # Create boundary condition
        self._conditions = uw.conditions.DirichletCondition(
            variable=self.TField,
            indexSetsPerDof=(self.top + self.bottom,))

        # Create Eq System
        self._system = uw.systems.SteadyStateHeat(temperatureField = self.TField,
                                                  fn_diffusivity = 1.0,
                                                  conditions = self._conditions)

        self._solver = uw.systems.Solver(self._system)

    def _solve_sle(self):
        self._solver.solve()

    def _advect_surface(self, dt):

        if self.Model.top_wall.data.size > 0:
            # Extract top surface
            x = self.Model.mesh.data[self.top.data][:, 0]
            y = self.Model.mesh.data[self.top.data][:, 1]

            # Extract velocities from top
            vx = self.Model.velocityField.data[self.top.data][:,0]
            vy = self.Model.velocityField.data[self.top.data][:,1]

            # Advect top surface
            x2 = x + vx * nd(dt)
            y2 = y + vy * nd(dt)

            # Spline top surface
            f = interp1d(x2, y2, kind='cubic')

            self.TField.data[self.top.data, 0] = f(x)

    def _update_mesh(self):

        with self.Model.mesh.deform_mesh():
            # Last dimension is the vertical dimension
            self.Model.mesh.data[:, -1] = self.TField.data[:, 0]

    def solve(self, dt):

        # First we advect the surface
        self._advect_surface(dt)

        uw.barrier()

        # Then we solve the system of linear equation
        self._solve_sle()

        # Finally we update the mesh
        self._update_mesh()
