Dynamic rc settings
===================

You can dynamically change the default rc settings in a python script or
interactively from the python shell. All of the rc settings are stored
in a dictionary-like variable called ``UWGeodynamics.rcParams``, which
is global to the UWGeodynamics package. rcParams can be modified
directly, for example::

.. code:: python

       import UWGeodynamics as GEO
       GEO.rcParams['solver'] = "mumps"
       GEO.rcParams['penalty'] = 1e6

The ``UWGeodynamics.rcdefaults`` command will restore the standard
UWGeodynamics default settings.

There is some degree of validation when setting the values of rcParams,
see ``UWGeodynamics.rcsetup`` for details.

The ``uwgeodynamicsrc`` file
============================

UWGeodynamics uses ``uwgeodynamicsrc`` configuration files to customize
all kinds of properties, which we call ``rc settings`` or
``rc parameters``. For now, you can control the defaults of a limited
set of property in matplotlib UWGeodynamics looks for
``uwgeodynamicsrc`` in four locations, in the following order:

1. ``uwgeodynamicsrc`` in the current working directory, usually used
   for specific customizations that you do not want to apply elsewhere.

2. ``$UWGEODYNAMICSRC`` if it is a file, else
   ``$UWGEODYNAMICSRC/uwgeodynamicsrc``.

3. It next looks in a user-specific place, depending on your platform:

   -  On Linux, it looks in ``.config/uwgeodynamics/uwgeodynamicsrc``
      (or ``$XDG_CONFIG_HOME/uwgeodynamics/uwgeodynamicsrc``) if you’ve
      customized your environment.

   -  On other platforms, it looks in
      ``.uwgeodynamics/uwgeodynamicsrc``.

4. ``{INSTALL}/UWGeodynamics/uwgeo-data/uwgeodynamicsrc``, where
   ``{INSTALL}`` is something like ``/usr/lib/python2.5/site-packages``
   on Linux, and maybe ``C:\\Python25\\Lib\\site-packages`` on Windows.
   Every time you install matplotlib, this file will be overwritten, so
   if you want your customizations to be saved, please move this file to
   your user-specific matplotlib directory.

To display where the currently active ``uwgeodynamicsrc`` file was
loaded from, one can do the following:

.. code:: python

     >>> import UWGeodynamics as GEO
     >>> GEO.uwgeodynamics_fname()
     '/home/foo/.config/uwgeodynamics/uwgeodynamicsrc'

See below for a sample.

\_uwgeodynamicsrc-sample:
-------------------------

A sample uwgeodynamicsrc file
-----------------------------
