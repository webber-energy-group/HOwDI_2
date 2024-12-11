# HOwDI (Hydrogen Optimization with Distribution Infrastructure)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

* This README is under development. For temporary reference, please consult [Modeling Hydrogen Infrastructure with the HOwDI Model](http://dx.doi.org/10.26153/tsw/43878)

## Installation

1. Create a conda environment from the `env.yml` file and activate it. This may be faster with [mamba](https://mamba.readthedocs.io/en/latest/).

    ```bash
    conda env create -f env.yml
    ```

    Note: M1 Mac users should change the Python version to 3.9 in `env.yml`

2. Activate the conda environment

    ```bash
    conda activate HOwDI
    ```

3. Install an editable version of HOwDI in your HOwDI environment with pip.

    ```bash
    pip install -e .
    ```

If you have already installed HOwDI but would like to update the dependencies:

```bash
conda env update --file env.yml --prune
```

## Usage

Within a directory that contains a subdirectory named "inputs" (that contains the necessary inputs), run the model:

```bash
(HOwDI) ~ ls
inputs
(HOwDI) ~ HOwDI run
```

Use `HOwDI run -h` for a list of options.

## Postprocessing Tools

HOwDI has several postprocessing tools. Use `HOwDI help` for a full list.

```bash
Create a figure:        HOwDI create_fig
Traceback:              HOwDI traceback
Traceforward:           HOwDI traceforward
```

## Contributing

HOwDI uses the Black code style. Please format your code accordingly before making a pull request.

## : Local Config

Adjust the config without worrying about git tracking by creating a file called `HOwDI/config_local.yml`. Add key `db` and follow with db path.


## Usage Continued

### Generating Hub Data

To create a new scenario and run the model from scratch, users will first have to generate a `hubs.csv` file that has columns called hub, type, latitude, longitude, and status. For the hub column, users should use camel case. If multiple hubs are located in the same city, users should use numbers to differentiate them. If the hub is located at a Love's station, users can add a 'L' to the end of the hub name. All hub names must be unique. 

For type, users can specify whether the hub is a city or a Love's station.

For status, 1 indicates a major hub, 0 a regular hub, and -1 a minor hub. Major hubs are more likely to connect to other hubs than regular or minor hubs.

In the same directory, users should create an `arcs_whitelist.csv` and an `arcs_blacklist.csv`. It is important that the hubs included in these files follow the same naming conventions as in `hubs.csv` so the model can match them up. 

Once you have created your `hubs.csv`, `arcs_whitelist.csv`, and `arcs_blacklist.csv` you are ready to create your hub data. An example of this command is: 
```bash
HOwDI create_hub_data -d hub_dir
```

`-d` signals that you will be providing a path to the `hub_dir` where `hub_dir` must contain `hubs.csv`, `arcs_whitelist.csv`, and `arcs_blacklist.csv`. Replace `hub_dir` with your path to your equivalent directory with the required files. From that command, `hubs.geojson`, `roads.csv`, and `arcs.csv` will be generated. `roads.csv` should include road geometry for every arc created in `arcs.csv`.

To create a figure of the generated arcs, use the following command:
```bash
HOwDI create_hub_data -f -shp US_COUNTY_SHPFILE/US_county_cont.shp -d hub_dir
```
`-d` again signals that you are providing a path to `hub_dir`. Replace `hub_dir` with your equivalent directory path with the input files. `-f` prompts the command to generate a figure. `-shp` signals that you will be providing a path to the shapefile you would like the figure to use. `US_COUNTY_SHPFILE/...` is an example of this path, but provide it with your shapefile path. 

To add the price multiplier values to `hubs.csv` according to their county, there is an automatic command. An example of this command is: 
```bash
HOwDI create_hub_data -pm ./cities_and_loces_pm.csv -d ./hub_dir -r -i ./hub_dir
```
`-pm` signals that you want to add price multiplier values to `hubs.csv` and the following path is your intended price multiplier file. Replace this path with the path to your price multiplier file. `-d` again signals that you are providing a path to `hub_dir`. Replace `hub_dir` with your equivalent directory path with the input files. `-r` signals that you desire to replace the `hubs.csv` and `arcs.csv` to use the new generated hubs and price multipliers while preserving other columns. `i` and the following directory provides the location of model inputs where files will be adjusted with new hubs. Replace this path with your desired location, although it can be helpful to use the same `hub_dir` directory or your equivalent. 

Use `HOwDI create_hub_data -h` for more detailed help descriptions.

### Creating an Input Scenario

Create an `inputs` directory for your scenario that includes:
```bash
arcs.csv
ccs.csv
conversion.csv
demand.csv
production_electric.csv
production_existing.csv
production_thermal.csv
settings.yml
distribution.csv
storage.csv
hubs.csv
```
Also within the `inputs` directory, create a `data` directory that includes:
```bash
hubs.geojson
roads.csv
```

In the `settings.yml` file, under '# optional options' set `hubs_dir` to the path of this `data` directory. If the `data` directory is within `inputs`, just the word "data" will suffice. 

Before running the model, make sure that hubs.csv includes the columns:
```bash
hub
type
latitude
longitude
status
build_smr-noCCS
build_smr-low
build_smr-high
build_smr-high-ptc
build_electrolyzer
build_electrolyzer-RE
transportationFuel_tonnesperday
transportationFuel_carbonSensitive_tonnesperday
industrialFuel_tonnesperday
existing_tonnesperday
ng_usd_per_mmbtu
e_usd_per_kwh
capital_pm
```

### Running Your Scenario

Within a directory that contains a subdirectory named "inputs" (that contains the necessary inputs), run the model:

```bash
(HOwDI) ~ ls
inputs
(HOwDI) ~ HOwDI run
```

Use `HOwDI run -h` for a list of options. Your new scenario should now have its corresponding outputs.