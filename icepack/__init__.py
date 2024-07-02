import glob
import yaml
import xarray as xr

def load_dataset(name, variable, config_path=None):
    """
    Loads and returns a .nc dataset into memory using the information in config.yaml

    Args:
        Name (str):         The name of the dataset in config.yaml
        Variable (str):     The name of the variable of interest (e.g., 'siconc', 'psl')
        config_path (str):  The path to config.yaml. If None, then assumed to be in pwd
        
    Returns:
        An xarray dataset, or a list of xarray datasets, depending on how many files there are for the dataset.
    """
    if type(config_path) == type(None):
        config_path = 'config.yaml'
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    paths = sorted(glob.glob(config['datasets'][name][variable]['path']))
    dataset = xr.open_dataset(paths[0]) if len(paths) == 1 else [xr.open_dataset(path) for path in paths]

    return dataset
