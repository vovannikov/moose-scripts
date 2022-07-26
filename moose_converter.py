import csv
import os, shutil
import time
import math
import sys
import numpy as np
import argparse
from csv import DictReader

def compute_initial_divisions(bottom_left, top_right, interface_width, elements_per_interface):

    domain_size   = top_right - bottom_left
    domain_width  = domain_size[0]
    domain_height = domain_size[1]
    domain_depth  = domain_size[2]

    initial_ny = 10
    initial_nx = int(domain_width / domain_height * initial_ny)
    initial_nz = int(domain_depth / domain_height * initial_ny)

    n_refinements = np.round(np.log2(elements_per_interface / interface_width *
                           domain_height / initial_ny))

    subdivisions = [initial_nx, initial_ny, initial_nz]
    
    return [subdivisions, n_refinements]

def append_line(input_data, line):

    input_data += line
    input_data += "\n"

    return input_data

def generate_moose_input(cloud_file, tpl_file, dim = 2, boundary_factor = 0.5, max_h_level_incr = 1, \
    t_end = 200., \
    A = 16., B = 1., L = 1., kappa_op = 0.5, kappa_c = 1.0, \
    Dvol = 1e-2, Dvap = 1e-10, Dsurf = 4.0, Dgb = 0.4):

    max_dim = 3

    input_data = ''
    try:
        with open(tpl_file, 'r') as tpl_file_handler:
            input_data = tpl_file_handler.read()
    except EnvironmentError: # parent of IOError, OSError *and* WindowsError where available
        print('Template file not found')

    particles = []

    with open(cloud_file, 'r') as read_obj:
        csv_dict_reader = DictReader(read_obj)
        for row in csv_dict_reader:

            p = {
                'c': [float(row['x']), float(row['y']), float(row['z'])],
                'r': float(row['r'])
            }
                
            particles.append(p)
    
    # Compute domain boundaries
    bottom_left = np.zeros(max_dim)
    top_right = np.zeros(max_dim)

    for d in range(dim):
        bottom_left[d] = sys.float_info.max
        top_right[d] = -sys.float_info.max

    r_max = -sys.float_info.max

    for p in particles:

        r_max = max(r_max, p['r'])

        for d in range(dim):
            if p['c'][d] - p['r'] < bottom_left[d]:
                bottom_left[d] = p['c'][d] - p['r']

            if p['c'][d] + p['r'] > top_right[d]:
                top_right[d] = p['c'][d] + p['r']

    for d in range(dim):
        bottom_left[d] -= boundary_factor * r_max
        top_right[d] += boundary_factor * r_max

    n_particles = len(particles)

    interface_width = 2.0
    elements_per_interface = 4
    subdivisions, n_refinements = compute_initial_divisions(bottom_left, top_right, interface_width, elements_per_interface)

    # Fill up the input file
    input_data = input_data.replace('{op_num}', str(n_particles))
    input_data = input_data.replace('{dim}', str(dim))

    for d in range(max_dim):
        placehoder_min = "min_{}".format(d)
        placehoder_max = "max_{}".format(d)

        input_data = input_data.replace('{' + placehoder_min + '}', "{:.4f}".format(bottom_left[d]))
        input_data = input_data.replace('{' + placehoder_max + '}', "{:.4f}".format(top_right[d]))

    input_data = input_data.replace('{nx}', str(subdivisions[0]))
    input_data = input_data.replace('{ny}', str(subdivisions[1]))
    input_data = input_data.replace('{nz}', str(subdivisions[2]))
    input_data = input_data.replace('{h_level}', str(n_refinements))
    input_data = input_data.replace('{max_h_level}', str(n_refinements + max_h_level_incr))

    var_name_base = 'gr'
    input_data = input_data.replace('{var_name_base}', str(var_name_base))
    op_list = [var_name_base + str(i) for i in range(n_particles)]
    gr_args_list = " ".join(op_list)
    op_coupling_list = ",".join(op_list)
    gr_kappa_list = " ".join(["kappa_op"]*n_particles)
    input_data = input_data.replace('{gr_args_list}', gr_args_list)
    input_data = input_data.replace('{gr_kappa_list}', gr_kappa_list)
    input_data = input_data.replace('{op_coupling_list}', op_coupling_list)

    # Free energy parameters
    input_data = input_data.replace('{energy_A}', str(A))
    input_data = input_data.replace('{energy_B}', str(B))
    input_data = input_data.replace('{energy_L}', str(L))
    input_data = input_data.replace('{energy_kappa_op}', str(kappa_op))
    input_data = input_data.replace('{energy_kappa_c}', str(kappa_c))

    # Mobility data
    input_data = input_data.replace('{Dvol}', str(Dvol))
    input_data = input_data.replace('{Dvap}', str(Dvap))
    input_data = input_data.replace('{Dsurf}', str(Dsurf))
    input_data = input_data.replace('{Dgb}', str(Dgb))

    # Simulation time
    input_data = input_data.replace('{t_end}', str(t_end))

    # Initial conditions section
    input_data += "[ICs]\n"

    x_positions = []
    y_positions = []
    z_positions = []
    radii = []

    counter = 0
    for p in particles:
        x_positions.append("{}".format(p['c'][0]))
        y_positions.append("{}".format(p['c'][1]))
        z_positions.append("{}".format(p['c'][2]))
        radii.append("{}".format(p['r']))

        input_data += "  [./ic_gr{}]\n".format(counter)
        input_data += "    int_width = {}\n".format(interface_width)
        input_data += "    x1 = {}\n".format(p['c'][0])
        input_data += "    y1 = {}\n".format(p['c'][1])
        input_data += "    z1 = {}\n".format(p['c'][2])
        input_data += "    radius = {}\n".format(p['r'])
        input_data += "    outvalue = 0.0\n"
        input_data += "    invalue = 1.0\n"
        input_data += "    variable = gr{}\n".format(counter)
        input_data += "    type = SmoothCircleIC\n"
        input_data += "  [../]\n"

        counter += 1

    input_data += "  [./multip]\n"
    input_data += "    x_positions = '{}'\n".format(" ".join(x_positions))
    input_data += "    y_positions = '{}'\n".format(" ".join(y_positions))
    input_data += "    z_positions = '{}'\n".format(" ".join(z_positions))
    input_data += "    radii = '{}'\n".format(" ".join(radii))
    input_data += "    int_width = {}\n".format(interface_width)
    input_data += "    3D_spheres = {}\n".format('true' if dim == max_dim else 'false')
    input_data += "    outvalue = 0.0\n"
    input_data += "    invalue = 1.0\n"
    input_data += "    variable = c\n"
    input_data += "    type = SpecifiedSmoothCircleIC\n"
    input_data += "  [../]\n"


    input_data += "[]\n"

    return input_data


# Main stuff

parser = argparse.ArgumentParser(description='Input options')
parser.add_argument("-c", "--cloud_file", type=str, help="Cloud file to output", required=True)
parser.add_argument("-s", "--save_file", type=str, help="Save file name", required=True)
parser.add_argument("-t", "--tpl_file", type=str, help="Template file", required=True)

args = parser.parse_args()

input_data = generate_moose_input(args.cloud_file, args.tpl_file)

with open(args.save_file, 'w') as f:
    f.write(input_data)

