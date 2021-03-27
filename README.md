# Seismic strength of the structure
### Description
This repository contains the structural calculations to determine the seismic strength of the building. The strength calclulations depend on the location of the building, soil conditions, and the dimensions of the building. The soil consitions and the location are related to probabilistic ground motion.
The ground motion coefficients are obtained from the USGS seismic website by sending a request to: https://earthquake.usgs.gov/designmaps). 
The tool is interactive and require multiple user inputs.
The code was developed in 2018 to automate calculations for a variety of different buildings.
The program is limited to a conventional (rectangular-shaped) buildings.

### Dependencies
 - python
 - numpy
 - urllib.parse
 - requests
 - matplotlib and matplotlib.pyplot

### Running Locally

* Install [Python](https://www.python.org/downloads/)


To run your application locally, run this inside the virtualenv:

```bash
python earthquake.py start
```

Your application will be running in your CLI.

# Contact details<a name="contact"></a>
If you have any questions, contact me via email: 

<a href="mailto:kathy.gomozova@gmail.com?"><img src="https://img.shields.io/badge/gmail-%23DD0031.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>