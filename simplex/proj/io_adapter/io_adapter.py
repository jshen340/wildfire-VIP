# io_adapter.py
import torch
import numpy as np

import os
from struct import pack

# all functions take 3D input, add the 3rd dim to 2D cell_counts and coordinates before writing

def create_folder(dirname):
    if (os.path.exists(dirname) == False):
        os.makedirs(dirname)


def Write_Grid(file_name, cell_counts, dx, domain_min):
    output = open(file_name, 'wb')
    try:
        output.write(cell_counts.tobytes())
        output.write(dx.tobytes())
        output.write(domain_min.tobytes())

    finally:
        output.close()


def Write_Particles(file_name, x, n, v=None, f=None, m=None):
    output = open(file_name, 'wb')
    try:
        output.write(pack('i', n))
        if (n != 0):
            placeholder_T = np.zeros(n * 3, dtype=np.float64)
            placeholder_T2 = np.zeros(n, dtype=np.float64)
            placeholder_i = np.zeros((n), dtype=int)
            output.write(x)
            if (v is not None):
                output.write(v)
            else:
                output.write(placeholder_T)
            if (f is not None):
                output.write(f)
            else:
                output.write(placeholder_T)
            if (m is not None):
                output.write(m)
            else:
                output.write(placeholder_T2)

            output.write(placeholder_T2)
            output.write(placeholder_i)

    finally:
        output.close()


def Write_Scalar_Field(file_name, s, counts):  # TODO
    output = open(file_name, 'wb')
    try:
        output.write(counts)
        output.write(s)
    finally:
        output.close()


def Write_Vector_Field(file_name, v, counts):
    output = open(file_name, 'wb')
    try:
        output.write(counts)
        output.write(v)
    finally:
        output.close()

def read():  # Compare
    output_dir = "D:/simplex/bin/win/io_adapter/Release\output"
    for i in range(1):
        frame_dir = os.path.join(output_dir, str(i))
        file_name = frame_dir + "/velocity"
        with open(file_name, 'rb') as f:
            a = f.read()
        # print(a)
        # file_name=frame_dir+"/particles"

        # file_name=frame_dir+"/phi"

        # file_name=frame_dir+"/velocity"


# read()

def test():
    output_dir = "output"
    for i in range(100):
        os.makedirs(output_dir+"/"+str(i))
    # grid property
    cell_counts = np.array([16, 16, 16], dtype=int)
    dx = 1 / cell_counts[0]
    domain_min = np.array([0, 0, 0], dtype=np.float64)

    # particle
    px = np.array([0, 0, 0, .5, .5, .5], dtype=np.float64)
    pn = 2
    vn = cell_counts[0] * cell_counts[1] * cell_counts[2]
    # scalar field property
    s = np.array(list(range(vn)), dtype=np.float64) / vn

    # vector field property
    v = np.zeros(vn * 3)
    for i in range(vn):
        v[i * 3] = .01
    # IO test
    for frame in range(100):
        px[3] += .1
        for i in range(vn):
            v[3 * i + 1] += .002

        frame_dir = os.path.join(output_dir, str(frame))

        file_name = frame_dir + "/grid"
        Write_Grid(file_name, cell_counts, dx, domain_min)

        file_name = frame_dir + "/particles"
        Write_Particles(file_name, px, pn)

        file_name = frame_dir + "/phi"
        # Write_Scalar_Field(file_name,s,cell_counts)

        file_name = frame_dir + "/velocity"
        Write_Vector_Field(file_name, v, cell_counts)

# test()