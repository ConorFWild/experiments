import yaml

ResidueID = tuple[str,str]
AtomID = tuple[str,str,str]

def read_yaml(path):
    with open(path, 'r') as f:
        dic = yaml.safe_load(f)

    return dic