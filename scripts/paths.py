import yaml
import os
import pathlib

def get_data_paths(mapping_path: str = '../../data_mapping.yml'):
    with open(mapping_path, 'r') as f:
        data_mapping = yaml.safe_load(f)

    root_dir = pathlib.Path(os.path.dirname(mapping_path))
    paths = {
        k: root_dir / data_mapping['folders']['data'] / p
        for k,p in data_mapping['datasets'].items()
    }

    return paths