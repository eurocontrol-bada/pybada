pyBADA documentation
==================================

pyBADA provides aircraft performance modelling, trajectory prediction and optimisation, and visualisation with BADA in Python.

To get started:

.. code-block:: bash

    pip install pybada


About BADA
==========

Efficient air traffic operations rely on the capability of air traffic management (ATM) systems to accurately predict aircraft trajectories. Likewise, ATM research and development activities require modelling and simulation tools capable of replicating real-life operations and aircraft performances as perfectly as possible.

To enable modelling and simulation of aircraft trajectories, every aircraft in operation shall have a corresponding Aircraft Performance Model (APM).

This is why EUROCONTROL, in cooperation with aircraft manufacturers and operating airlines, has created and maintains an APM called BADA (Base of Aircraft Data).

Owning to this cooperation, BADA prevails as a unique product provided by EUROCONTROL to the aviation community worldwide

BADA APM is based on a kinetic approach to aircraft performance modelling. It is made of two components: the Model specifications, which provide the theoretical fundaments used to calculate aircraft performance parameters, and the Datasets containing the aircraft-specific coefficients necessary to perform calculations. 



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pyBADA
   architecture
   userManual
   history


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`