"""
Author: Zhinuo J. Wang 5 July 2016

Use principal component analysis (PCA) to extract geometric modes and modal coefficients for a population of left
ventricular (LV) models segmented from patient MRI. These modes and modal coefficients will then be used in a framework
for estimating the reference geometry (i.e. estimating the modal coefficients which combine to re-create the LV
geometry) as well as the passive stiffness parameter simultaneously.

Input(s):
    - 28 patient LV models in 16 tricubic Hermite elements.
Output(s):
    - All geometric modes which describe variation in geometry within the population.
    - Corresponding eigenvalues/modal coefficients.
    - Small number of modes capable of describing 99% of the geometric variation.

Code for conversion of tricubic Hermite to tricubic Lagrange adapted from MATLAB code from Xiani Yan.
"""

import numpy as np
import os
np.set_printoptions(threshold=np.nan)

def _read_ipnode(filename):
    nodes = []
    # Read ipnode file. self.node has indices meaning: node_num, component, version, derivative_num
    with open(filename, 'r') as f:
        data = f.readlines()
    node_versions = []
    i = 12
    while i <= len(data):
        i += 1
        v = int(data[i].split()[-1])
        node_versions.append(v)
        comp = []
        i += 1
        for c in [1, 2, 3]:
            version = []
            for j in range(0, v):
                if v > 1:
                    i += 1
                derivative = []
                for count in [1, 2, 3, 4, 5, 6, 7, 8]:
                    derivative.append(float(data[i].split()[-1]))
                    i += 1
                version.append(derivative)
            i += 1
            comp.append(version)
        nodes.append(comp)
    return nodes

def _read_exelem(filename):
    elems = []
    scale_factors = np.zeros((40, 8))
    # Read elements and scale factors from .exelem file.
    with open(filename, 'r') as f:
        data = f.readlines()
    # Special elements 1 to 4.
    lines = [[563, 576], [670, 683], [777, 790], [884, 897]]
    for k in range(0, 4):
        temp = [int(j) for j in data[lines[k][0] - 2].split()]
        if k == 3:
            elems.append([34 + k, temp[0], temp[1], temp[2], 37 + k, temp[3], temp[4], temp[5]])
        elif k == 0:
            elems.append([temp[0], 35, temp[1], temp[2], temp[3], 38, temp[4], temp[5]])
        else:
            elems.append([35 + k - 1, 35 + k, temp[1], temp[2], 38 + k - 1, 38 + k, temp[4], temp[5]])
        temp = []
        for i in range(lines[k][0], lines[k][1]):
            temp.append([float(j) for j in data[i].split()])
        temp = sum(temp, [])
        for i in range(0, 8):
            scale_factors[elems[-1][i] - 1] = temp[i * 8:i * 8 + 8]
    # All other elements.
    i = 980
    while i < len(data):
        i += 9
        elems.append([int(j) for j in data[i].split()])
        i += 2
        temp = []
        for k in range(0, 13):
            temp.append([float(j) for j in data[k + i].split()])
        i += 12
        temp = sum(temp, [])
        for l in range(0, 8):
            scale_factors[elems[-1][l] - 1] = temp[l * 8:l * 8 + 8]
        i += 1
    return elems, scale_factors

def ipnode_to_fem_format(nodes, elements, scale_factors):
    print 'Converting from ipnode format to global DOFs vectors'
    #num_nodes_per_elem = 8
    num_derivs = 8
    num_nodes = len(nodes)
    #num_elems = len(elements)
    #num_params_per_elem = num_derivs * num_nodes_per_elem

    # Save nodal DOFs read in from ipnode file into global DOF vector - add on 6 duplicated nodes to cater for
    # collapsed nodes at apex.
    global_dofs_map = np.zeros((num_nodes+6, num_derivs)) # Mapping from node number and derivative format
    # to global number of DOFs format.
    global_dofs_param = np.zeros(((num_nodes+6)*num_derivs, 3)) # All global DOFs in x, y and z.
    n_global_dof = 0
    for n in range(0, num_nodes):
        for d in range(0, num_derivs):
            global_dofs_map[n, d] = n_global_dof
            for c in range(0, 3):
                global_dofs_param[n_global_dof, c] = nodes[n][c][0][d] #* scale_factors[n][d]
            n_global_dof += 1

    # Additional nodes for node 1 at apex:
    n_global_dof_end = num_nodes*num_derivs
    for i in range(1, 4):
        for d in range(0, num_derivs):
            global_dofs_map[i+33, d] = n_global_dof_end
            for c in range(0, 3):
                global_dofs_param[n_global_dof_end, c] = nodes[0][c][i][d] #* scale_factors[i + 33][d]
            n_global_dof_end +=1
    # Additional nodes for node 18 at apex:
    for i in range(1, 4):
        for d in range(0, num_derivs):
            global_dofs_map[i + 36, d] = n_global_dof_end
            for c in range(0, 3):
                global_dofs_param[n_global_dof_end, c] = nodes[17][c][i][d] #* scale_factors[i + 36][d]
            n_global_dof_end += 1

    """
    # Connectivity for Hermite elements - latent capability.
    # Global DOFs to elements mapping.
    global_dofs_local_elem_map = np.zeros((num_nodes_per_elem*num_derivs, num_elems))
    # Elements 1 to 4 are special apical elements that require node numbering using additional nodes.
    apical_elements = [[1, 35, 2, 3, 18, 38, 19, 20],
                       [35, 36, 3, 4, 38, 39, 20, 21],
                       [36, 37, 4, 5, 39, 40, 21, 22],
                       [37, 1, 5, 2, 40, 18, 22, 19]]
    for e in range(0, 4):
        n_local_dof = 0
        for n in range(0, num_nodes_per_elem):
            for d in range(0, num_derivs):
                global_node_number = apical_elements[e][n]
                global_dofs_local_elem_map[n_local_dof, e] = global_dofs_map[global_node_number-1, d]
                n_local_dof += 1
    # All other elements are in the format read in from the exelem file.
    for e in range(4, num_elems):
        n_local_dof = 0
        for n in range(0, num_nodes_per_elem):
            for d in range(0, num_derivs):
                global_node_number = elements[e][n]
                global_dofs_local_elem_map[n_local_dof, e] = global_dofs_map[global_node_number-1, d]
                n_local_dof += 1

    # Global dofs to global element map.
    global_dofs_global_elem_map = np.zeros((num_params_per_elem*num_elems, n_global_dof))
    global_elem_params = np.zeros((num_params_per_elem * num_elems, 3))
    row_number = 0  # Rows count through the total DOFs in each element consecutively. The columns count through the
    # total number of DOFs of the whole mesh.
    for e in range(0, num_elems):
        for n_local_dof in range(0, num_params_per_elem):
            global_dof_number = global_dofs_local_elem_map[n_local_dof, e]
            global_dofs_global_elem_map[row_number, global_dof_number] = 1
            global_elem_params[row_number, 0]= global_dofs_param[global_dof_number, 0]
            global_elem_params[row_number, 1] = global_dofs_param[global_dof_number, 1]
            global_elem_params[row_number, 2] = global_dofs_param[global_dof_number, 2]
            row_number += 1
    """
    return global_dofs_param

def fem_format_to_ipnode(global_dofs_param):
    print 'Convert from global DOFs vectors into ipnode format.'
    nodes = []
    for n in range(0, 34):
        if n == 0:
            component = []
            for c in range(0, 3):
                version = []
                for v in range(0, 4):
                    deriv = []
                    for d in range(0, 8):
                        if v < 1:
                            deriv.append(global_dofs_param[d, c])
                        else:
                            deriv.append(global_dofs_param[(v + 33)*8 + d, c])
                    version.append(deriv)
                component.append(version)
        elif n == 17:
            component = []
            for c in range(0, 3):
                version = []
                for v in range(0, 4):
                    deriv = []
                    for d in range(0, 8):
                        if v < 1:
                            deriv.append(global_dofs_param[17*8 + d, c])
                        else:
                            deriv.append(global_dofs_param[(v + 36)*8 + d, c])
                    version.append(deriv)
                component.append(version)
        else:
            component = []
            for c in range(0, 3):
                version = []
                deriv = []
                for d in range(0, 8):
                    deriv.append(global_dofs_param[n*8+d, c])
                version.append(deriv)
                component.append(version)
        nodes.append(component)
    return nodes

def hermite_to_bezier(fem_hermite, h_to_b_matrix):
    print 'Converting Hermite DOFs to Bezier DOFs.'
    fem_bezier = np.zeros(fem_hermite.shape) # Initialise Bezier parameters matrix using shape of Hermite matrix.
    num_derivs = 8
    # Convert the Hermite DOFs to Bezier DOFs - being consistent with ordering of derivatives in the global_dofs_param.
    for c in range(0, 3):
        n_dof = 0
        for n in range(0, 40):
            fem_hermite_temp = fem_hermite[n_dof:n_dof + 8, c]
            for d_row in range(0, num_derivs):
                temp = 0
                for d_column in range(0, num_derivs):
                    temp += fem_hermite_temp[d_column] * h_to_b_matrix[d_row, d_column]
                fem_bezier[n_dof, c] = temp
                n_dof += 1
    return fem_bezier

def bezier_to_hermite(fem_bezier, h_to_b_matrix):
    print 'Converting Bezier DOFs to Hermite DOFs'
    fem_hermite = np.zeros(fem_bezier.shape)  # Initialise Hermite parameters matrix using shape of Bezier matrix.
    num_derivs = 8
    Hg = np.linalg.inv(h_to_b_matrix)  # Inverse the h-to-b conversion matrix to convert from B to H.
    # Convert the Bezier DOFs to Hermite DOFs - being consistent with ordering of derivatives in the global_dofs_param.
    for c in range(0, 3):
        n_dof = 0
        for n in range(0, 40):
            fem_bezier_temp = fem_bezier[n_dof:n_dof + 8, c]
            for d_row in range(0, num_derivs):
                temp = 0
                for d_column in range(0, num_derivs):
                    temp += fem_bezier_temp[d_column] * Hg[d_row, d_column]
                fem_hermite[n_dof, c] = temp
                n_dof += 1
    return fem_hermite

def export_geometry(nodes, filename):
    print 'Exporting current geometry to '+filename + '.ipnode and '+filename + '.exnode files.'
    # Write out ipnode file.
    with open(filename + '.ipnode', 'w') as f:
        f.write(' CMISS Version 2.1  ipnode File Version 2\n')
        f.write(' Heading: Elements created in Perl\n\n')
        f.write(' The number of nodes is [    34]:     34\n')
        f.write(' Number of coordinates [3]: 3\n')
        f.write(' Do you want prompting for different versions of nj=1 [N]? Y\n')
        f.write(' Do you want prompting for different versions of nj=2 [N]? Y\n')
        f.write(' Do you want prompting for different versions of nj=3 [N]? Y\n')
        f.write(' The number of derivatives for coordinate 1 is [0]: 7\n')
        f.write(' The number of derivatives for coordinate 2 is [0]: 7\n')
        f.write(' The number of derivatives for coordinate 3 is [0]: 7\n\n')

        for n in range(0, len(nodes)):
            f.write(' Node number [     ' + str(n + 1) + ']:      ' + str(n + 1) + '\n')
            if n in [0, 17]:
                for c in range(0, 3):
                    f.write(' The number of versions for nj='+str(c+1)+' is [1]:  4\n')
                    for v in range(0, 4):
                        f.write(' For version number '+str(v+1)+':\n')
                        f.write(' The Xj(' + str(c + 1) + ') coordinate is [ 0.66029E+02]:    '
                                + str(nodes[n][c][v][0]) + '\n')
                        f.write(' The derivative wrt direction 1 is [ 0.00000E+00]:  '
                                + str(nodes[n][c][v][1]) + '\n')
                        f.write(' The derivative wrt direction 2 is [-0.74612E-01]:  '
                                + str(nodes[n][c][v][2]) + '\n')
                        f.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:  '
                                + str(nodes[n][c][v][3]) + '\n')
                        f.write(' The derivative wrt direction 3 is [ 0.10000E+01]:   '
                                + str(nodes[n][c][v][4]) + '\n')
                        f.write(' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:    '
                                + str(nodes[n][c][v][5]) + '\n')
                        f.write(' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:    '
                                + str(nodes[n][c][v][6]) + '\n')
                        f.write(' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:    '
                                + str(nodes[n][c][v][7]) + '\n')
                f.write('\n')
            else:
                for c in range(0, 3):
                    f.write(' The number of versions for nj='+str(c+1)+' is [1]:  1\n')
                    f.write(' The Xj('+str(c+1)+') coordinate is [ 0.66029E+02]:    '+str(nodes[n][c][0][0])+'\n')
                    f.write(' The derivative wrt direction 1 is [ 0.00000E+00]:  '+str(nodes[n][c][0][1])+'\n')
                    f.write(' The derivative wrt direction 2 is [-0.74612E-01]:  '+str(nodes[n][c][0][2])+'\n')
                    f.write(' The derivative wrt directions 1 & 2 is [ 0.00000E+00]:  '+str(nodes[n][c][0][3])+'\n')
                    f.write(' The derivative wrt direction 3 is [ 0.10000E+01]:   '+str(nodes[n][c][0][4])+'\n')
                    f.write(' The derivative wrt directions 1 & 3 is [ 0.00000E+00]:    '+str(nodes[n][c][0][5])+'\n')
                    f.write(' The derivative wrt directions 2 & 3 is [ 0.00000E+00]:    '+str(nodes[n][c][0][6])+'\n')
                    f.write(' The derivative wrt directions 1, 2 & 3 is [ 0.00000E+00]:    '+str(nodes[n][c][0][7])+'\n')
                f.write('\n')
    # Export to exnode file.
    with open(filename + '.exnode', 'w') as f:
        f.write(' Group name: ReferenceWallModel\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 4\n')
        f.write('   y.  Value index=33, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 4\n')
        f.write('   z.  Value index=65, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 4\n')
        f.write(' Node:            1\n')
        for c in range(0, 3):
            for v in range(0, 4):
                for d in range(0, 8):
                    f.write('\t'+str(nodes[0][c][v][d]))
            f.write('\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)\n')
        f.write('   y.  Value index= 9, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)\n')
        f.write('   z.  Value index=17, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)\n')
        for n in range(1, 17):
            f.write(' Node:            '+str(n+1)+'\n')
            for c in range(0, 3):
                for d in range(0, 8):
                    f.write('\t'+str(nodes[n][c][0][d]))
                f.write('\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 4\n')
        f.write('   y.  Value index=33, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 4\n')
        f.write('   z.  Value index=65, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3), #Versions= 4\n')
        f.write(' Node:            18\n')
        for c in range(0, 3):
            for v in range(0, 4):
                for d in range(0, 8):
                    f.write('\t' + str(nodes[17][c][v][d]))
            f.write('\n')
        f.write(' #Fields=1\n')
        f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
        f.write('   x.  Value index= 1, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)\n')
        f.write('   y.  Value index= 9, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)\n')
        f.write('   z.  Value index=17, #Derivatives= 7 (d/ds1,d/ds2,d2/ds1ds2,d/ds3,d2/ds1ds3,d2/ds2ds3,d3/ds1ds2ds3)\n')
        for n in range(18, len(nodes)):
            f.write(' Node:            ' + str(n + 1) + '\n')
            for c in range(0, 3):
                for d in range(0, 8):
                    f.write('\t' + str(nodes[n][c][0][d]))
                f.write('\n')
    # Export to exelem file - 1.0 scale factors since they are now included in the nodes.

def main():
    # Import tricubic Hermite model from CIM models.
    node_file = 'DSWall.ipnode' # Temporarily
    # hard-coded.
    elem_file = 'DSWall.exelem'
    nodes = _read_ipnode(node_file)
    [elements, scale_factors] = _read_exelem(elem_file)

    # Read in conversion matrices
    h_to_b_matrix = np.loadtxt('h_to_b_matrix.txt', dtype='float')
    #BH = np.loadtxt('BH.txt', dtype='float') # For mapping connectivity from Bezier to Hermite. - latent capability
    #HB = np.loadtxt('HB.txt', dtype='float') # For mapping connectivity from Hermite to Bezier. - latent capability

    # Convert into finite element format for tricubic Hermite
    fem_hermite = ipnode_to_fem_format(nodes, elements, scale_factors)

    # Convert tricubic Hermite parameters into tricubic Bezier parameters.
    fem_bezier = hermite_to_bezier(fem_hermite, h_to_b_matrix)

    # Just to check, convert Bezier back to Hermite and output to exnode file.
    fem_hermite_back = bezier_to_hermite(fem_bezier, h_to_b_matrix)
    nodes_back = fem_format_to_ipnode(fem_hermite_back)
    export_geometry(nodes_back, 'DSWall_back')

    print 'Completed without crashing.'

if __name__ == '__main__':
    main()