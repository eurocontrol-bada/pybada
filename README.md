# pyBADA

<a href="https://github.com/eurocontrol/pybada/blob/main/LICENCE.txt"><img alt="License: EUPL" src="https://img.shields.io/badge/license-EUPL-3785D1.svg"></a>
![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB.svg?logo=python&logoColor=white)
<a href="https://github.com/eurocontrol/pybada"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>

The BADA aircraft performance toolbox for Python

To get started

```bash
pip install pyBADA
```

## Examples

-   `file_parser`: example of BADA file parser and retrieval of some basic BADA parameters for all BADA3/4/H
-   `optimum_speed_altitude`: example of calculation of optimum speeds and altitude for BADA4 and BADAH aircraft
-   `ac_trajectory`: a simple, but complete aircraft trajectory for BADA3 and BADA4 aircraft
-   `ac_trajectory_GPS`: an example of a simple, but complete aircraft trajectory for BADA3 and BADA4 aircraft including geodesic calculations

## Development

```bash
# Optionally, set up a virtual env and activate it
python3 -m venv env
source env/bin/activate
# Install package in editable mode
pip install -e .
# Install a couple of packages for formatting and linting
pip install -r requirements-dev.txt
```

## License

BADA and pyBADA are developed and maintained by [EUROCONTROL](https://www.eurocontrol.int/).

This project is licensed under the European Union Public License v1.2 - see the [LICENSE](https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12) file for details.
