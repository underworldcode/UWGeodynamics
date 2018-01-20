#!/usr/bin/python
import h5py
import os
import numpy as np
import xml.etree.ElementTree as ET
import argparse

__author__ = "Romain Beucher"


def convert_to_ascii(file_input, file_output=None):

    root_dir = os.path.dirname(file_input)

    with open(file_input, "r") as f:
        tree = ET.parse(f)
        root = tree.getroot()

    grid = root[0].find("Grid")
    geometry = grid.find("Geometry")
    attribute = grid.find("Attribute")

    file_geometry = geometry[0].text
    file_attribute = attribute[0].text

    f = h5py.File(os.path.join(root_dir, file_geometry.split(":")[0]), "r")
    geometry = f["/data"]

    f = h5py.File(os.path.join(root_dir, file_attribute.split(":")[0]), "r")
    data = f["/data"]

    new = np.zeros((geometry.value.shape[0], geometry.value.shape[1] + 1))
    new[:, :-1] = geometry.value
    new[:, -1] = data.value.ravel()

    if not file_output:
        file_output = os.path.split(file_input)[-1].split(".")[0] + ".asc"

    np.savetxt(file_output, new)


def main():

    description = """Convert UW .h5 file to ASCII
                   for import into MOVE/PETREL"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    args = parser.parse_args()

    file_input = str(args.input)
    file_output = str(args.output)

    convert_to_ascii(file_input, file_output)
