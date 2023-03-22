import os
import json


class DotDict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def _read(path):
    with open(path) as json_file:
        obj = DotDict(json.load(json_file))
    return obj


def read_titles():
    return _read(os.path.dirname(__file__) + "/../config/titles.json")


def read_erros():
    return _read(os.path.dirname(__file__) + "/../config/errors.json")
