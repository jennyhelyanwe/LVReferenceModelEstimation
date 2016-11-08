import numpy as np
import os


# Helper functions for performing the PCA on the study cohort.
# Zhinuo Wang
# Procrustean transform code adapted from code by Mahyar Osanlouy 1 Nov 2016

def read_ipnode(filename):
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


def export_geometry(nodes, filename):
    #print 'Exporting current geometry to '+filename + '.ipnode and '+filename + '.exnode files.'
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


def ipnode_to_dofs_format(nodes):
    num_derivs = 8
    num_nodes = len(nodes)
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
                global_dofs_param[n_global_dof, c] = nodes[n][c][0][d]
            n_global_dof += 1
    # Additional nodes for node 1 at apex:
    n_global_dof_end = num_nodes*num_derivs
    for i in range(1, 4):
        for d in range(0, num_derivs):
            global_dofs_map[i+33, d] = n_global_dof_end
            for c in range(0, 3):
                global_dofs_param[n_global_dof_end, c] = nodes[0][c][i][d]
            n_global_dof_end +=1
    # Additional nodes for node 18 at apex:
    for i in range(1, 4):
        for d in range(0, num_derivs):
            global_dofs_map[i + 36, d] = n_global_dof_end
            for c in range(0, 3):
                global_dofs_param[n_global_dof_end, c] = nodes[17][c][i][d]
            n_global_dof_end += 1
    return global_dofs_param


def dofs_to_ipnode_format(global_dofs_param):
    #print 'Convert from global DOFs vectors into ipnode format.'
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


def convert_arithmetic_to_unit_sf(input_filename, output_filename):
    with open('ConvertToUnitSF.com', 'r') as f:
        data = f.readlines()
    data[5] = 'fem def node;r;'+input_filename+'\n'
    data[6] = 'fem def elem;r;'+input_filename+'\n'
    data[11] = 'fem exp node;'+output_filename+'\n'
    data[12] = 'fem exp elem;'+output_filename+'\n'
    data[13] = 'fem def node;w;'+output_filename+'\n'
    data[14] = 'fem def elem;w;'+output_filename+'\n'
    with open('ConvertToUnitSF.com', 'w') as f:
        f.writelines(data)
    os.system('cm ConvertToUnitSF.com')


def hermite_dofs_to_bezier_dofs(fem_hermite, h_to_b_matrix):
    # print 'Converting Hermite DOFs to Bezier DOFs.'
    fem_bezier = np.zeros(fem_hermite.shape)  # Initialise Bezier parameters matrix using shape of Hermite matrix.
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


def bezier_dofs_to_hermite_dofs(fem_bezier, h_to_b_matrix):
    #print 'Converting Bezier DOFs to Hermite DOFs'
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


def procrustes(reference_model, model):
    X = reference_model.reshape([320, 3])
    Y = model.reshape([320, 3])

    n, m = np.shape(X)
    ny, my = np.shape(Y)

    # Evaluate average of reference and input DOFs.
    muX = np.average(X, axis=0)
    muY = np.average(Y, axis=0)

    # Offset DOFs by average value.
    X0 = X - muX
    Y0 = Y - muY

    # Evaluate sum of squares for offset DOFs
    ssX = (X0**2.).sum()
    ssY = (Y0**2.).sum()

    # Evaluate centred Frobenius norm
    normX = np.sqrt(ssX)
    normY = np.sqrt(ssY)

    # Scale to unit norms
    X0 /= normX
    Y0 /= normY

    if my < m:
        Y0 = np.concatenate((Y0, np.zeros(n, m-my)),0) # padding input unit norm DOFs

    A = np.dot(X0.T, Y0)
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    V = Vt.T
    T = np.dot(V, U.T)

    traceTA = s.sum()

    # Find optimum scaling of Y
    b = traceTA * normX / normY

    # Standardised distance between X and b*Y*T + c
    d = 1 - traceTA**2

    # Transform input DOFs.
    Z = normX*traceTA*np.dot(Y0, T) + muX

    # Evaluate translation matrix
    c = muX - b*np.dot(muY, T)
    tform = {'rotation':T, 'scale':b, 'translation':c}
    return d, Z, tform


def procrustes_reverse(procrusted_model, tform):
    nodes = procrusted_model.reshape([320, 3])
    translated_node = nodes - np.array(tform['translation'])
    #scaled_node = translated_node / tform['scale'] Do not scale for better PCA.
    rotated_node = np.dot(translated_node,np.linalg.inv(tform['rotation']))
    return rotated_node


def evaluate_reconstruction_error(study_name):
    with open('CompareReconstructedModel.com', 'r') as f:
        data = f.readlines()
    data[6] = 'fem def node;r;Models_UnitSF/DSWall_'+study_name+'\n'
    data[7] = 'fem def elem;r;Models_UnitSF/DSWall_'+study_name+'\n'

    data[12] = 'fem def data;w;Models_UnitSF/DSWall_'+study_name+'_SP_Endo\n'
    data[13] = 'fem exp data;Models_UnitSF/DSWall_'+study_name+'_SP_Endo\n'
    data[18] = 'fem def data;w;Models_UnitSF/DSWall_' + study_name + '_SP_Epi\n'
    data[19] = 'fem exp data;Models_UnitSF/DSWall_' + study_name + '_SP_Epi\n'

    data[29] = 'fem def node;r;ReconstructionModels/Reconstructed_'+study_name+'\n'
    data[30] = 'fem def elem;r;Models_UnitSF/DSWall_'+study_name+'\n'

    data[32] = 'fem def data;r;Models_UnitSF/DSWall_'+study_name+'_SP_Endo\n'
    data[35] = 'fem list data;ReconstructionModels/'+study_name+'_Endo_Error error\n'
    data[36] = 'fem export data;ReconstructionModels/'+study_name+'_Endo_Error as EndoError error\n'

    data[38] = 'fem def data;r;Models_UnitSF/DSWall_' + study_name + '_SP_Epi\n'
    data[41] = 'fem list data;ReconstructionModels/' + study_name + '_Epi_Error error\n'
    data[42] = 'fem export data;ReconstructionModels/' + study_name + '_Epi_Error as EpiError error\n'
    with open('CompareReconstructedModel.com', 'w') as f:
        f.writelines(data)
    os.system('cm CompareReconstructedModel.com -> ComLog.out')
