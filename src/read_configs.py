import json

class DotDict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def _read(path):
    with open(path) as json_file:
        obj = DotDict(json.load(json_file))
    return obj

def read_titles(path='./src/config/titles.json'):
    return _read(path)    

def read_erros(path= './src/config/errors.json'):
    return _read(path)