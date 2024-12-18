Getting started
==================================

.. code-block:: bash

    pip install pyBADA


Running on unsupported environments
-----------------------------------

You won't receive support for it, but you can pass the flag ``--ignore-requires-python`` to install pyBADA on an unsupported Python version.


Providing custom models
-----------------------

To load your own BADA models, when instantiating an aircraft, pass as ``filePath`` the parent folder containing the desired BADA models. For example:

E.g. If your ``/home/<USER>/bada-models/BADA4/4.3/`` folder contains:

.. code-block:: text

    A320-212
    .
    ├── A320-212.ATF
    ├── A320-212_final_.mat
    ├── A320-212_ISA+20.PTD
    ├── A320-212_ISA+20.PTF
    ├── A320-212_ISA.PTD
    ├── A320-212_ISA.PTF
    ├── A320-212.xml
    ├── ECON.OPT
    ├── LRC.OPT
    ├── MEC.OPT
    ├── MRC.OPT
    └── OPTALT.OPT
    A320-232
    .
    ├── A320-232.ATF
    ├── A320-232_final_.mat
    ├── A320-232_ISA+20.PTD
    ├── A320-232_ISA+20.PTF
    ├── A320-232_ISA.PTD
    ├── A320-232_ISA.PTF
    ├── A320-232.xml
    ├── ECON.OPT
    ├── LRC.OPT
    ├── MEC.OPT
    ├── MRC.OPT
    └── OPTALT.OPT

You can instantiate an A320-232 like this:

.. code-block:: python

    AC = Bada4Aircraft(
        badaVersion=badaVersion,
        acName="A320-232",
        filePath="/home/<USER>/bada-models/BADA4/4.3",
    )

In the case of BADA3, this means directly pointing to the folder containing the model files.


