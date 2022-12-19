from configparser import ConfigParser
from pathlib import Path


def read_ini(path: str):
    config = ConfigParser()
    config.read(path)

    # Create same ini without comments
    ini_raw_path = (Path.cwd() / '_tmp' / 'raw_config.ini').resolve()
    config_raw = ConfigParser()
    for section in config.sections():
        config_raw[section] = {}
        for (key, val) in config.items(section):
            config_raw[section][key] = supress_comment_from_value(val)
    with open(ini_raw_path, 'w') as configfile:
        config_raw.write(configfile)

    return config_raw


def supress_comment_from_value(value: str) -> str:
    return value.split('#')[0].rstrip()
