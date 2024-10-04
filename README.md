# pyBADA

<a href="https://github.com/eurocontrol/pybada/blob/main/LICENCE.txt"><img alt="License: EUPL" src="https://img.shields.io/badge/license-EUPL-3785D1.svg"></a>
<a href="https://pypi.org/project/pyBADA"><img alt="Released on PyPi" src="https://img.shields.io/pypi/v/pyBADA.svg"></a>
![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB.svg?logo=python&logoColor=white)
<a href="https://github.com/eurocontrol/pybada"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>

The BADA aircraft performance toolbox for Python

To get started

```bash
pip install pyBADA
```

## Examples

-   `file_parser`: BADA file parser and retrieval of some basic BADA parameters for all BADA3/4/H
-   `BADAData`: loading complete BADA3 dataset and retrievel of a specific parameters for a specific aircraft
-   `optimum_speed_altitude`: calculation of optimum speeds and altitude for BADA4 and BADAH aircraft
-   `ac_trajectory`: simple, but complete aircraft trajectory for BADA3 and BADA4 aircraft
-   `ac_trajectory_GPS`: simple, but complete aircraft trajectory for BADA3 and BADA4 aircraft including geodesic calculations
-   `heli_trajectory`: simple, but complete helicopter trajectory for BADAH aircraft

## Development

```bash
# Optionally, set up a virtual env and activate it
python3 -m venv env
source env/bin/activate
# Install package in editable mode
pip install -e .
# Install a couple of packages for formatting, linting and building the docs
pip install -e .[dev]
# To build the docs
cd docs
make html
```

## License

BADA and pyBADA are developed and maintained by [EUROCONTROL](https://www.eurocontrol.int/).

This project is licensed under the European Union Public License v1.2 - see the [LICENSE](https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12) file for details.
