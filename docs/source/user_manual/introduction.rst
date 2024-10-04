Introduction
============

Identification
--------------
PyBADA is a Python-based implementation of BADA APM for aircraft performance modelling, trajectory prediction and optimisation with BADA, including methods to calculate any of the aircraft’s performance characteristics (lift, drag, thrust, and fuel flow), operational flight envelope, BADA airline procedure model and calculation of PTF and PTD files, as defined by the BADA User Manual. It includes modules to calculate atmospheric properties, unit conversions (including the speed conversion) and basic geodesic calculations.

PyBADA offers a full implementation of BADA3, BADA4 and BADAH. 

The complete implementation then builds on top of BADA performance calculation by providing several trajectory segment calculations for all phases of flight (cruise, climb and descend) including acceleration and deceleration computation.

EUROCONTROL is the owner of pyBADA

Purpose
-------
PyBADA library has been created with several use cases in mind:

- Help users starting with BADA and aircraft performance modelling to better understand how such a model can be implemented
- Lower the workload of users planning to implement BADA in their tools, by reusing an already exiting fully functional library;
- Minimizing the errors and unnecessary mistakes while implementing BADA performance calculations;
- Help users who are building their own tools to validate the results against the library implemented and validated by the BADA team at EUROCONTROL;
- Provide a starting platform for research organisations and universities, to build new and interesting tools utilizing BADA performance modelling, without requiring too much time to be spent on BADA implementation, and rather focusing on the research part.

Glossary of Acronyms
--------------------

.. list-table:: Glossary of Acronyms
   :widths: 15 85
   :header-rows: 1

   * - Acronym
     - Explanation
   * - ARPM
     - Airline Procedure Model
   * - BADA
     - Base of Aircraft Data
   * - BADA3
     - Base of Aircraft Data Family 3
   * - BADA4
     - Base of Aircraft Data Family 4
   * - BADAH
     - Base of Aircraft Data Family H
   * - PTD
     - Performance Table Data
   * - PTF
     - Performance Table File
   * - pyBADA
     - Python Implementation of BADA
   * - SW
     - Software
   * - TAS
     - True Air Speed
   * - TBP
     - Turboprop
   * - TCL
     - Trajectory Computation Light