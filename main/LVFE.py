import sys
import os
from opencmiss.iron import iron
import numpy
from numpy import array
import exfile
from util_surface_data import readIpdata
from bkcolours import *

#sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'], 'cm', 'bindings', 'python')))

# This script solves finite elasticity models of the left ventricle of the heart.
# Author: ZJW 2015-2016

class LV:
    def __init__(self):
        #self.compNodes = iron.ComputationalNumberOfNodesGet()
        #self.compNodeNum = iron.ComputationalNodeNumberGet()
        self.compNodes = 1
        self.compNodeNum = 1
        #os.system('pwd')
        #iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL, [1, 2, 3, 4, 5], 'diagnostic', ["FiniteElasticity_GaussDeformationGradientTensor"])
        self.numOfXi = 3
        self.nodes = []
        self.elems = []
        self.fibre = []
        self.scaleFactors = numpy.zeros((40, 8))
        self.scalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN


    def setup(self):
        # Coordiante
        print 'Coordinates and regions and basis setup'
        self._setup_coordinates()
        self._setup_regions()
        self._setup_bases()

        # Geometric set up:
        print 'Geometric setup'
        self._read_ref_model()
        self._setup_cavity_model()

        # Wall reference geometry
        print 'wall mesh, decomposition, geometry setup'
        self._setup_wall_mesh()
        self._setup_wall_decomposition()
        self._setup_wall_geometry()
        self._setup_fibre()
        # Cavity reference geometry
        print 'cavity mesh, decomposition, geometry setup'
        self._setup_cavity_mesh()
        self._setup_cavity_decomposition()
        self._setup_cavity_geometry()

        # Mechanics set up.
        self.endoPressure = 0.0
        self.TCa = 0.0
        self.endoIncrem = 0.0
        self.activeIncrem = 0.0
        print 'Define mechanics'
        self._define_mechanics()
        self._setup_equations_set()
        self._setup_material()
        self._setup_dependent()
        self._setup_equations()
        self._setup_deformed()
        self._setup_hydro()
        self.activeCurrent = 0.0
        self.endoCurrent = 0.0
        self.endoNodes = [1, 35, 36, 37, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        self.epiNodes = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34]
        self.epiBaseNodes = [31, 32, 33, 34]
        self.baseBC = array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        self._export_wall("ref")
        self._export_cavity("ref")

    def solve(self, warm_solve_toggle):
        # Export reference models first.
        self._update_cavity()
        #print self.endoPressureLoadSteps
        for i in range(1, len(self.endoPressureLoadSteps)):
            self.endoPressure = self.endoPressureLoadSteps[i]
            self.baseBC = array(self.baseBCSteps[i])
            print bkcolours.OKBLUE + 'Set endocardial pressure to: ' + str(self.endoPressure) + ' kPa' + bkcolours.ENDC
            print bkcolours.OKBLUE + 'Set epicardial basal nodes coordinates to:' + bkcolours.ENDC
            print self.baseBC
            if self.TCa > 0:
                self._update_independent()
            self._setup_problem()
            self._setup_boundary_conditions(warm_solve_toggle)
            self.problem.Solve()
            self._update_cavity()
            self.problem.Finalise()
            self.solverEquations.Finalise()
            self._update_cavity()
        self._export_wall("def")
        self._export_cavity("def")

    def finalise(self):
        self.cavityCoordinateSystem.Finalise()
        self.wallCoordinateSystem.Finalise()
        self.cavityRegion.Finalise()
        self.wallRegion.Finalise()
        self.cavityBasis.Finalise()
        self.geomBasis.Finalise()
        self.pressureBasis.Finalise()

    def def_material(self, c1):
        stiffness = [c1, 8.60937380179457, 3.66592870916253, 25.7671807349653]
        for c, param in enumerate(stiffness, 1):
            self.materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                           iron.FieldParameterSetTypes.VALUES, c, param)
        print bkcolours.HEADER + 'Current C1: '+str(c1) + bkcolours.ENDC

    def initialise_dependent(self, current_frame):
        self._read_def_model(current_frame)
        # Create reference geometry for nodes which don't need versions.
        for n in range(0, len(self.hydroDOFs)):
            for c in range(0, 3):
                    for d in range(0, 8):
                        self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES, 1, d+1,
                                                                     n+1, c+1, self.deformedDOFs[n][c][d])

        temp = numpy.concatenate(self.deformedScaleFactors)
        self.dependentField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.U, 1, temp)
        self.fibreField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.U, 1, temp)
        # Set hydrostatic pressure
        for n in range(1, len(self.hydroDOFs)):
            self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                         iron.FieldParameterSetTypes.VALUES, 1, 1, n, 4,
                                                         self.hydroDOFs[n-1])
        self.dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
        self.dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

    def initialise_dependent_fromlegacy(self, ipinit_filename, elem_filename, node_filename):
        for i in range(0, len(self.deformedDOFs)):
            for k in range(0, 3):
                    for j in range(0, 8):
                        self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES, 1, j+1,
                                                                     i+1, k+1, self.deformedDOFs[i][k][0][j])
        # Initialise reference positions for the other three versions of node 1.
        #print self.deformedDOFs[0][1]
        for v in range(1, 4):
            for c in range(0, 3):
                for d in range(0, 8):
                        self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES, 1, d+1,
                                                                     v+34, c+1, self.deformedDOFs[0][c][v][d])
        # Initialise reference positions for the other three versions of node 18.
        for i in range(1, 4):
            for k in range(0, 3):
                for j in range(0, 8):
                    for v in range(1, 4):
                        self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES, 1, j+1,
                                                                     i+37, k+1, self.deformedDOFs[17][k][i][j])
        temp = numpy.concatenate(self.deformedScaleFactors)
        self.dependentField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.U, 1, temp)
        self.fibreField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.U, 1, temp)
        # Set hydrostatic pressure
        for n in range(1, len(self.hydroDOFs)+1):
            self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                         iron.FieldParameterSetTypes.VALUES, 1, 1, n, 4,
                                                         self.hydroDOFs[n-1])
        for n in range(35, 38):
            self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                         iron.FieldParameterSetTypes.VALUES, 1, 1, n, 4,
                                                         self.hydroDOFs[0])
        for n in range(38, 41):
            self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                         iron.FieldParameterSetTypes.VALUES, 1, 1, n, 4,
                                                         self.hydroDOFs[17])
        self.dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
        self.dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

    def projection(self, filename_epi, filename_endo, frame_num):
        print bkcolours.OKGREEN + 'Projecting data to surface...' + bkcolours.ENDC
        # Epi error.
        print 'Epicardial projection...'
        epi_data = readIpdata(filename_epi)
        epi_data_points = iron.DataPoints()
        epi_data_points.CreateStart(self.wallRegion, len(epi_data))
        for n in range(0, int(len(epi_data))):
            epi_data_points.ValuesSet(n+1, epi_data[n][:])
        epi_data_points.CreateFinish()
        epiDataProjectionUserNumber = 1
        epi_data_projection = iron.DataProjection()
        epi_data_projection.CreateStart(epiDataProjectionUserNumber, epi_data_points, self.wallMesh)
        epi_data_projection.ProjectionTypeSet(iron.DataProjectionProjectionTypes.BOUNDARY_FACES)
        epi_data_projection.ProjectionCandidatesSet([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
                                                    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6])
        epi_data_projection.CreateFinish()
        epi_data_projection.DataPointsProjectionEvaluate(self.dependentField)
        epi_error_vector = []
        for n in range(0, len(epi_data)):
            # Reversing sign of projection vector due to bug in source code.
            epi_error_vector.append(-epi_data_projection.ResultProjectionVectorGet(n+1, 3))
        epi_SE = []
        for n in range(0, len(epi_data)):
            epi_SE.append(epi_data_projection.ResultDistanceGet(n+1)**2)  # Squared error Euclidean distance.
        epi_sum_SE = sum(epi_SE)
        epi_MSE = epi_sum_SE/len(epi_data)
        epi_RMSE = numpy.sqrt(epi_MSE)
        print 'Epi MSE of frame ' + str(frame_num) + ' from OpenCMISS results is ' + str(epi_MSE)
        print 'Epi RMSE of frame ' + str(frame_num) + ' from OpenCMISS results is ' + str(epi_RMSE)
        self._export_error('fitting_error/EpiProjectionError_'+str(frame_num)+'.exdata', epi_data_points,
                           epi_error_vector)
        epi_data_points.Destroy()

        # Endo error.
        print 'Endocardial projection...'
        endo_data = readIpdata(filename_endo)
        endo_data_points = iron.DataPoints()
        endo_data_points.CreateStart(self.wallRegion, len(endo_data))
        for n in range(0, int(len(endo_data))):
            endo_data_points.ValuesSet(n+1, endo_data[n][:])
        endo_data_points.CreateFinish()
        endoDataProjectionUserNumber = 2
        endo_data_projection = iron.DataProjection()
        endo_data_projection.CreateStart(endoDataProjectionUserNumber, endo_data_points, self.wallMesh)
        endo_data_projection.ProjectionTypeSet(iron.DataProjectionProjectionTypes.BOUNDARY_FACES)
        endo_data_projection.ProjectionCandidatesSet([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
                                                     [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3])
        endo_data_projection.CreateFinish()
        endo_data_projection.DataPointsProjectionEvaluate(self.dependentField)
        endo_error_vector = []
        for n in range(0, len(endo_data)):
            # Reversing sign of projection vector due to bug in source code.
            endo_error_vector.append(-endo_data_projection.ResultProjectionVectorGet(n+1, 3))

        endo_SE = []
        for n in range(0, len(endo_data)):
            endo_SE.append(endo_data_projection.ResultDistanceGet(n+1)**2)  # Squared error Euclidean distance.
        endo_sum_SE = sum(endo_SE)   # Sum of squared error.
        endo_MSE = endo_sum_SE/len(endo_data) # Mean of sum of squared error: mean squared error.
        endo_RMSE = numpy.sqrt(endo_MSE)
        print 'Endo MSE of frame ' + str(frame_num) + ' from OpenCMISS results is ' + str(endo_MSE)
        print 'Endo RMSE of frame ' + str(frame_num) + ' from OpenCMISS results is ' + str(endo_RMSE)
        self._export_error('fitting_error/EndoProjectionError_'+str(frame_num)+'.exdata', endo_data_points,
                           endo_error_vector)
        endo_data_points.Destroy()

        mse = (epi_sum_SE + endo_sum_SE)/(len(epi_data) + len(endo_data))   # Mean squared error for both surfaces
        rmse = numpy.sqrt(mse)
        print 'Current total MSE for frame '+str(frame_num)+' is '+str(mse)
        print 'Current total RMSE for frame '+str(frame_num)+' is '+str(rmse)
        return mse

    def set_pressure(self, pressure, pressure_prev):
        self.endoIncrem = pressure - pressure_prev
        endoIncremLimit = 0.2
        self.loadsteps = numpy.ceil(abs(self.endoIncrem / endoIncremLimit))
        steps = numpy.linspace(0, self.endoIncrem, self.loadsteps+1)
        self.endoPressureLoadSteps = []
        for i in range(0, len(steps)):
            self.endoPressureLoadSteps.append(pressure_prev + steps[i])

        if self.endoIncrem == 0:  # Means we are doing a warm solve
            self.endoPressureLoadSteps = [0, pressure]

    def set_base_displacement(self, nodeDataIndices, cur_epi, next_epi, cur_endo, next_endo, cur_frame_num):
        nodes = self._get_node_displacement(nodeDataIndices, cur_epi, next_epi, cur_endo, next_endo)
        nodes = array(nodes)
        base_coords = self._get_base_coords('ForwardSolveSolution/LVWall_'+str(cur_frame_num)+'.exnode')
        x_disp = self._get_ave_x_displacement(nodeDataIndices, next_epi, base_coords)
        numSteps = self.loadsteps + 1
        #print 'Start base coordinates'
        #print base_coords
        temp = base_coords + nodes[4:8][:]
        for n in range(0, 4):
            temp[n][0] = base_coords[n][0] + x_disp
        #print 'Final basal coordinates should be:'
        #print temp
        self.baseBCSteps = []
        self.baseBCSteps.append(base_coords.tolist())
        increm = nodes[4:8][:]/(numSteps -1)
        x_increm = x_disp/(numSteps - 1)
        for i in range(0, int(numSteps-1)):
            #print i
            #print self.baseBCSteps[i]
            epi_nodes = (self.baseBCSteps[i] + increm).tolist()
            for n in range(0, 4):
                epi_nodes[n][0] = self.baseBCSteps[i][n][0] + x_increm
            self.baseBCSteps.append(epi_nodes)
            #print 'Epi nodes at loadstep ' + str(i + 1)
            #print epi_nodes
        #epi_nodes = base_coords + nodes[4:8][:]

        #for n in range(0, 4):
        #    epi_nodes[n][0] = base_coords[n][0] + x_disp
        #self.baseBC = epi_nodes
        #print self.baseBC
        #print 'Take steps:'
        #print self.baseBCSteps[0]
        #print self.baseBCSteps[1]
        #print self.baseBCSteps[2]
        #quit()
        #print self.baseBCSteps
        #quit()

    def set_base_displacement_zero(self):
        # Set zero displacement at epicardial base.
        base_coords = self._get_base_coords('LVWall_ref.part0.exnode')
        self.baseBCSteps = []
        self.baseBCSteps.append(base_coords.tolist())
        numSteps = self.loadsteps + 1
        for i in range(0, int(numSteps -1)):
            self.baseBCSteps.append(base_coords.tolist())

    def set_t_ca(self, t_ca):
        self.TCa = t_ca
        self.activeIncremLimit = 3.0

    def _export_error(self, filename, data, vector):
        with open(filename, 'w') as f:
            f.write(' Group name: ProjectionError\n')
            f.write(' #Fields=2\n')
            f.write(' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
            f.write('   x.  Value index= 1, #Derivatives=0\n')
            f.write('   y.  Value index= 2, #Derivatives=0\n')
            f.write('   z.  Value index= 3, #Derivatives=0\n')
            f.write(' 2) error, field, rectangular cartesian, #Components=3\n')
            f.write('   x.  Value index= 4, #Derivatives=0\n')
            f.write('   y.  Value index= 5, #Derivatives=0\n')
            f.write('   z.  Value index= 6, #Derivatives=0\n')
            for n in range(0, len(vector)):
                temp = data.ValuesGet(n+1, 3)
                f.write(' Node:     '+str(n+1)+'\n')
                f.write('   '+str(temp[0])+'\t'+str(temp[1])+'\t'+str(temp[2]) +
                        '\n   '+str(vector[n][0])+'\t'+str(vector[n][1]) +
                        '\t'+str(vector[n][2])+'\n')

    def _read_def_model(self, current_frame):
        print bkcolours.OKBLUE + 'Reading in warm start solution at frame ' + str(current_frame) + bkcolours.ENDC
        self.warmNodes = exfile.Exnode("LVInflation_"+str(current_frame)+".exnode")
        self.warmElems = exfile.Exelem("LVInflation_"+str(current_frame)+".exelem")
        self.deformedScaleFactors = self.warmElems.scale_factors
        hydro = []
        for node_num in range(1, self.warmNodes.num_nodes + 1):
            value = self.warmNodes.node_value("HydrostaticPressure", "1", node_num, 1)
            hydro.append(value)
        self.hydroDOFs = hydro
        coords = []
        for n in range(1, self.warmNodes.num_nodes+1):
            temp = []
            for c in [1, 2, 3]:
                component_name = ["x", "y", "z"][c-1]
                deriv = []
                for d in [1, 2, 3, 4, 5, 6, 7, 8]:
                    value = self.warmNodes.node_value("DeformedGeometry", component_name, n, d)
                    deriv.append(value)
                temp.append(deriv)
            coords.append(temp)
        coords = tuple(coords)
        self.deformedDOFs = coords

    def _read_def_model_fromlegacy(self, ipinit_filename, elem_filename, node_filename):
        nodesX = []
        nodesY = []
        nodesZ  = []
        try:
            fid = open(ipinit_filename, 'r')
        except IOError:
            print 'ERROR: Unable to open ', ipinit_filename
            return
        try:
            junk = fid.readline()
            toggle = True
            while toggle:
                if junk == " Dependent variable initial conditions:\n":
                    toggle = False
                else:
                    junk = fid.readline()
            junk = fid.readline()
            junk = fid.readline()

            # Read DOFs in X direction
            temp = fid.readline()
            temp = temp.split()
            num = temp[len(temp)-1]
            while num != '0':
                if (num == '1') | (num == '18'):
                    versions = []
                    for v in [1, 2, 3, 4]:
                        junk= fid.readline()
                        derivs = []
                        for i in [0, 1, 2, 3, 4, 5, 6, 7]:
                            temp = fid.readline()
                            temp = temp.split()
                            #print temp
                            derivs.append(float(temp[len(temp)-1]))
                        versions.append(derivs)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                    nodesX.append(versions)
                else:
                    derivs = []
                    versions = []
                    for i in [0, 1, 2, 3, 4, 5, 6, 7]:
                        temp = fid.readline()
                        temp = temp.split()
                        #print temp
                        derivs.append(float(temp[len(temp)-1]))
                    versions.append(derivs)
                    nodesX.append(versions)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
            #nodesX = numpy.array(nodesX)

            # Read DOFs in Y direction
            junk = fid.readline()
            junk = fid.readline()
            temp = fid.readline()
            temp = temp.split()
            num = temp[len(temp)-1]
            while num != '0':
                if (num == '1') | (num == '18'):
                    versions = []
                    for v in [1, 2, 3, 4]:
                        junk= fid.readline()
                        derivs = []
                        for i in [0, 1, 2, 3, 4, 5, 6, 7]:
                            temp = fid.readline()
                            temp = temp.split()
                            derivs.append(float(temp[len(temp)-1]))
                        versions.append(derivs)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                    nodesY.append(versions)
                else:
                    derivs = []
                    versions = []
                    for i in [0,1,2,3,4,5,6,7]:
                        temp = fid.readline()
                        temp = temp.split()
                        derivs.append(float(temp[len(temp)-1]))
                    versions.append(derivs)
                    nodesY.append(versions)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
            #nodesY = numpy.array(nodesY)

            # Read DOFs in Y direction
            junk = fid.readline()
            junk = fid.readline()
            temp = fid.readline()
            temp = temp.split()
            num = temp[len(temp)-1]
            while num != '0':
                if (num == '1') | (num == '18'):
                    versions = []
                    for v in [1, 2, 3, 4]:
                        junk= fid.readline()
                        derivs = []
                        for i in [0, 1, 2, 3, 4, 5, 6, 7]:
                            temp = fid.readline()
                            temp = temp.split()
                            derivs.append(float(temp[len(temp)-1]))
                        versions.append(derivs)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
                    nodesZ.append(versions)
                else:
                    derivs = []
                    versions = []
                    for i in [0,1,2,3,4,5,6,7]:
                        temp = fid.readline()
                        temp = temp.split()
                        derivs.append(float(temp[len(temp)-1]))
                    versions.append(derivs)
                    nodesZ.append(versions)
                    temp = fid.readline()
                    temp = temp.split()
                    num = temp[len(temp)-1]
            #nodesZ = numpy.array(nodesZ)

            # The indices for nodes goes: component (x,y,z), node number, derivative number.
            nodes = []
            for i in range(0, len(nodesX)):
                comp = []
                comp.append(nodesX[i])
                comp.append(nodesY[i])
                comp.append(nodesZ[i])
                nodes.append(comp)

            # Read hydrostatic pressure at nodes
            junk = fid.readline()
            junk = fid.readline()

            node_idx = 0
            temp = fid.readline()
            temp = temp.split()
            num = temp[len(temp)-1]
            hydro = []
            while num!='0':
                temp = fid.readline()
                temp = temp.split()
                hydro.append(float(temp[len(temp)-1]))
                temp = fid.readline()
                temp = temp.split()
                num = temp[len(temp)-1]
            hydro = numpy.array(hydro)
            self.deformedDOFs = nodes
            self.hydroDOFs = hydro
        finally:
            fid.close()
        with open(elem_filename, 'r') as f:
            data = f.readlines()
        # Special elements 1 to 4.
        elems = []
        scaleFactors = numpy.zeros((40, 8))
        lines = [[670, 683], [886, 899], [1102, 1115], [1318, 1331]]
        for k in range(0, 4):
            temp = [int(j) for j in data[lines[k][0]-2].split()]

            if k == 3:
                elems.append([34 + k, temp[0], temp[1], temp[2], 37+k, temp[3],  temp[4], temp[5]])
            elif k == 0:
                elems.append([temp[0], 35, temp[1], temp[2], temp[3], 38, temp[4], temp[5]])
            else:
                elems.append([35+k-1, 35 + k, temp[1], temp[2], 38+k-1, 38 + k, temp[4],  temp[5]])
            temp = []
            for i in range(lines[k][0], lines[k][1]):
                temp.append([float(j) for j in data[i].split()])
            temp = sum(temp, [])
            for i in range(0, 8):
                scaleFactors[elems[-1][i]-1] = temp[i*8:i*8+8]
        # All other elements.
        i = 1523
        while i < len(data):
            i += 9
            elems.append([int(j) for j in data[i].split()])
            i += 2
            temp = []
            for k in range(0, 13):
                temp.append([float(j) for j in data[k+i].split()])
            i += 14
            temp = sum(temp, [])
            for l in range(0, 8):
                scaleFactors[elems[-1][l]-1] = temp[l*8:l*8+8]
            i += 1
        self.deformedScaleFactors = scaleFactors

    def _get_node_displacement(self, Node_Data, cur_epi, next_epi, cur_endo, next_endo):
        epi = readIpdata(cur_epi)
        epi_def = readIpdata(next_epi)
        # Calculate the nodal displacement during inflation
        Node31_x_def = epi_def[Node_Data[0], 0] - epi[Node_Data[0], 0]
        Node32_x_def = epi_def[Node_Data[1], 0] - epi[Node_Data[1], 0]
        Node33_x_def = epi_def[Node_Data[2], 0] - epi[Node_Data[2], 0]
        Node34_x_def = epi_def[Node_Data[3], 0] - epi[Node_Data[3], 0]

        Node31_y_def = epi_def[Node_Data[0], 1] - epi[Node_Data[0], 1]
        Node32_y_def = epi_def[Node_Data[1], 1] - epi[Node_Data[1], 1]
        Node33_y_def = epi_def[Node_Data[2], 1] - epi[Node_Data[2], 1]
        Node34_y_def = epi_def[Node_Data[3], 1] - epi[Node_Data[3], 1]

        Node31_z_def = epi_def[Node_Data[0], 2] - epi[Node_Data[0], 2]
        Node32_z_def = epi_def[Node_Data[1], 2] - epi[Node_Data[1], 2]
        Node33_z_def = epi_def[Node_Data[2], 2] - epi[Node_Data[2], 2]
        Node34_z_def = epi_def[Node_Data[3], 2] - epi[Node_Data[3], 2]

        print 'Nodal displacement: '
        print 'Node 31 x, y, z diplacements = ', Node31_x_def, Node31_y_def, Node31_z_def
        print 'Node 32 x, y, z diplacements = ', Node32_x_def, Node32_y_def, Node32_z_def
        print 'Node 33 x, y, z diplacements = ', Node33_x_def, Node33_y_def, Node33_z_def
        print 'Node 34 x, y, z diplacements = ', Node34_x_def, Node34_y_def, Node34_z_def

        Node31_def = [Node31_x_def, Node31_y_def, Node31_z_def]
        Node32_def = [Node32_x_def, Node32_y_def, Node32_z_def]
        Node33_def = [Node33_x_def, Node33_y_def, Node33_z_def]
        Node34_def = [Node34_x_def, Node34_y_def, Node34_z_def]

        # Read in the data at ED and ES
        endo = readIpdata(cur_endo)
        endo_def = readIpdata(next_endo)

        print '================ Write out ipinit file based on displacement  ===================='
        # Calculate the nodal displacement during inflation
        Node14_x_def = endo_def[Node_Data[0], 0] - endo[Node_Data[0], 0]
        Node15_x_def = endo_def[Node_Data[1], 0] - endo[Node_Data[1], 0]
        Node16_x_def = endo_def[Node_Data[2], 0] - endo[Node_Data[2], 0]
        Node17_x_def = endo_def[Node_Data[3], 0] - endo[Node_Data[3], 0]

        Node14_y_def = endo_def[Node_Data[0], 1] - endo[Node_Data[0], 1]
        Node15_y_def = endo_def[Node_Data[1], 1] - endo[Node_Data[1], 1]
        Node16_y_def = endo_def[Node_Data[2], 1] - endo[Node_Data[2], 1]
        Node17_y_def = endo_def[Node_Data[3], 1] - endo[Node_Data[3], 1]

        Node14_z_def = endo_def[Node_Data[0], 2] - endo[Node_Data[0], 2]
        Node15_z_def = endo_def[Node_Data[1], 2] - endo[Node_Data[1], 2]
        Node16_z_def = endo_def[Node_Data[2], 2] - endo[Node_Data[2], 2]
        Node17_z_def = endo_def[Node_Data[3], 2] - endo[Node_Data[3], 2]

        print 'Nodal displacement: '
        print 'Node 14 x, y, z diplacements = ', Node14_x_def, Node14_y_def, Node14_z_def
        print 'Node 15 x, y, z diplacements = ', Node15_x_def, Node15_y_def, Node15_z_def
        print 'Node 16 x, y, z diplacements = ', Node16_x_def, Node16_y_def, Node16_z_def
        print 'Node 17 x, y, z diplacements = ', Node17_x_def, Node17_y_def, Node17_z_def

        Node14_def = [Node14_x_def, Node14_y_def, Node14_z_def]
        Node15_def = [Node15_x_def, Node15_y_def, Node15_z_def]
        Node16_def = [Node16_x_def, Node16_y_def, Node16_z_def]
        Node17_def = [Node17_x_def, Node17_y_def, Node17_z_def]

        return Node14_def, Node15_def, Node16_def, Node17_def, Node31_def, Node32_def, Node33_def, Node34_def
        """
        node31 = epi_def[Node_Data[0], :]
        node32 = epi_def[Node_Data[1], :]
        node33 = epi_def[Node_Data[2], :]
        node34 = epi_def[Node_Data[3], :]
        return node31, node32, node33, node34
        """

    def _get_ave_x_displacement(self, Node_Data, next_epi, base_coords):
        CAPEpi_def = readIpdata(next_epi)
        # Calculate the nodal displacement during inflation
        Node31_x_def = CAPEpi_def[Node_Data[0], 0] - base_coords[0, 0] # x displacement measured against model.
        Node32_x_def = CAPEpi_def[Node_Data[1], 0] - base_coords[1, 0]
        Node33_x_def = CAPEpi_def[Node_Data[2], 0] - base_coords[2, 0]
        Node34_x_def = CAPEpi_def[Node_Data[3], 0] - base_coords[3, 0]
        x_def = array([Node31_x_def, Node32_x_def, Node33_x_def, Node34_x_def])
        mean_x = numpy.mean(x_def)
        print 'Mean x displacement: '+str(mean_x)
        return mean_x

    def _get_base_coords(self, filename):
        allnodes = exfile.Exnode(filename)
        basenodes = []
        for n in [31, 32, 33, 34]:
            coords = []
            for c in ["x", "y", "z"]:
                coords.append(allnodes.node_value("DeformedGeometry", c, n, 1))
            basenodes.append(coords)
        basenodes = array(basenodes)
        return basenodes

    def _read_ref_model(self):
        self.nodes = []
        self.elems = []
        self.fibre = []
        # Read ipnode file. self.node has indices meaning: node_num, component, version, derivative_num
        with open("ReferenceWall.ipnode", 'r') as f:
            data = f.readlines()
        self.node_versions = []
        i = 12
        while i <= len(data):
            i += 1
            v = int(data[i].split()[-1])
            self.node_versions.append(v)
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
            self.nodes.append(comp)
        # Read elements and scale factors from .exelem file.
        with open("ReferenceWall.exelem", 'r') as f:
            data = f.readlines()
        # Special elements 1 to 4.
        lines = [[563, 576], [670, 683], [777, 790], [884, 897]]
        for k in range(0, 4):
            temp = [int(j) for j in data[lines[k][0]-2].split()]
            if k == 3:
                self.elems.append([34 + k, temp[0], temp[1], temp[2], 37+k, temp[3],  temp[4], temp[5]])
            elif k == 0:
                self.elems.append([temp[0], 35, temp[1], temp[2], temp[3], 38, temp[4], temp[5]])
            else:
                self.elems.append([35+k-1, 35 + k, temp[1], temp[2], 38+k-1, 38 + k, temp[4],  temp[5]])
            temp = []
            for i in range(lines[k][0], lines[k][1]):
                temp.append([float(j) for j in data[i].split()])
            temp = sum(temp, [])
            for i in range(0, 8):
                self.scaleFactors[self.elems[-1][i]-1] = temp[i*8:i*8+8]
        # All other elements.
        i = 980
        while i < len(data):
            i += 9
            self.elems.append([int(j) for j in data[i].split()])
            i += 2
            temp = []
            for k in range(0, 13):
                temp.append([float(j) for j in data[k+i].split()])
            i += 12
            temp = sum(temp, [])
            for l in range(0, 8):
                self.scaleFactors[self.elems[-1][l]-1] = temp[l*8:l*8+8]
            i += 1
        # Read fibre orientation.
        with open('DTMRI_CIMModel.ipfibr', 'r') as f:
            data = f.readlines()
        i = 20
        while i <= len(data):
            i += 1
            comp = []
            for c in [1, 2, 3]:
                deriv = []
                for d in [1, 2, 3, 4, 5, 6, 7, 8]:
                    deriv.append(float(data[i].split()[-1]))
                    i += 1
                comp.append(deriv)
            self.fibre.append(comp)
            i += 1
        # Identify apex nodes.
        self.apexEndoNodes = [1, 35, 36, 37]
        self.apexEpiNodes = [18, 38, 39, 40]

    def _setup_cavity_model(self):
        apex = [self.nodes[0][0][0][0], self.nodes[0][1][0][0], self.nodes[0][2][0][0]]
        base_nodes_x = [self.nodes[13][0][0][0], self.nodes[14][0][0][0], self.nodes[15][0][0][0],
                        self.nodes[16][0][0][0]]
        base_nodes_y = [self.nodes[13][1][0][0], self.nodes[14][1][0][0], self.nodes[15][1][0][0],
                        self.nodes[16][1][0][0]]
        base_nodes_z = [self.nodes[13][2][0][0], self.nodes[14][2][0][0], self.nodes[15][2][0][0],
                        self.nodes[16][2][0][0]]
        self.baseMean = [numpy.mean(base_nodes_x), numpy.mean(base_nodes_y), numpy.mean(base_nodes_z)]
        base_apex = apex[0] - self.baseMean[0]

    def _setup_coordinates(self):
        # Set up coordinate system for the wall.
        self.wallCoordinateSystemUserNumber = 1
        self.wallCoordinateSystem = iron.CoordinateSystem()
        self.wallCoordinateSystem.CreateStart(self.wallCoordinateSystemUserNumber)
        self.wallCoordinateSystem.DimensionSet(self.numOfXi)
        self.wallCoordinateSystem.CreateFinish()
        # Set up coordinate system for the cavity.
        self.cavityCoordinateSystemUserNumber = 2
        self.cavityCoordinateSystem = iron.CoordinateSystem()
        self.cavityCoordinateSystem.CreateStart(self.cavityCoordinateSystemUserNumber)
        self.cavityCoordinateSystem.DimensionSet(self.numOfXi)
        self.cavityCoordinateSystem.CreateFinish()

    def _setup_regions(self):
        # Set up wall region.
        self.wallRegionUserNumber = 1
        self.wallRegion = iron.Region()
        self.wallRegion.CreateStart(self.wallRegionUserNumber, iron.WorldRegion)
        self.wallRegion.CoordinateSystemSet(self.wallCoordinateSystem)
        self.wallRegion.LabelSet("Wall")
        self.wallRegion.CreateFinish()
        # Set up LV cavity region.
        self.cavityRegionUserNumber = 2
        self.cavityRegion = iron.Region()
        self.cavityRegion.CreateStart(self.cavityRegionUserNumber, iron.WorldRegion)
        self.cavityRegion.CoordinateSystemSet(self.cavityCoordinateSystem)
        self.cavityRegion.LabelSet("Cavity")
        self.cavityRegion.CreateFinish()

    def _setup_bases(self):
        # Wall geometry basis (tricubic Hermite)
        self.geomBasisUserNumber = 1
        self.geomBasis = iron.Basis()
        self.geomBasis.CreateStart(self.geomBasisUserNumber)
        self.geomBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
        self.geomBasis.NumberOfXiSet(self.numOfXi)
        self.geomBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*self.numOfXi)
        self.geomBasis.QuadratureNumberOfGaussXiSet([4]*self.numOfXi)
        self.geomBasis.QuadratureLocalFaceGaussEvaluateSet(True)
        self.geomBasis.CreateFinish()
        # Pressure basis (trilinear)
        self.pressureBasisUserNumber = 2
        self.pressureBasis = iron.Basis()
        self.pressureBasis.CreateStart(self.pressureBasisUserNumber)
        self.pressureBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
        self.pressureBasis.NumberOfXiSet(self.numOfXi)
        self.pressureBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*self.numOfXi)
        self.pressureBasis.QuadratureNumberOfGaussXiSet([4]*self.numOfXi)
        self.pressureBasis.QuadratureLocalFaceGaussEvaluateSet(True)
        self.pressureBasis.CreateFinish()
        # Cavity basis (tricubic Hermite)
        self.cavityBasisUserNumber = 3
        self.cavityBasis = iron.Basis()
        self.cavityBasis.CreateStart(self.cavityBasisUserNumber)
        self.cavityBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
        self.cavityBasis.NumberOfXiSet(self.numOfXi)
        self.cavityBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*self.numOfXi)
        self.cavityBasis.QuadratureNumberOfGaussXiSet([4]*self.numOfXi)
        self.cavityBasis.QuadratureLocalFaceGaussEvaluateSet(True)
        self.cavityBasis.CreateFinish()

    def _setup_wall_mesh(self):
        # Set up wall mesh.
        self.wallMeshUserNumber = 1
        self.wallMesh = iron.Mesh()
        self.wallMesh.CreateStart(self.wallMeshUserNumber, self.wallRegion, self.numOfXi)
        self.wallMesh.NumberOfComponentsSet(2)
        self.wallMesh.NumberOfElementsSet(16)
        self.wallNodes = iron.Nodes()
        # Using 4 nodes at endo-cardial and epicardial apex, we need a total of 34 + 6 nodes.
        self.wallNodes.CreateStart(self.wallRegion, 40)
        self.wallNodes.CreateFinish()
        # Set up 16 element mesh for wall geometry.
        self.cubicWallElementsMeshComponentNumber = 1
        self.cubicWallElements = iron.MeshElements()
        self.cubicWallElements.CreateStart(self.wallMesh, self.cubicWallElementsMeshComponentNumber, self.geomBasis)
        for e in range(4, len(self.elems)):
            self.cubicWallElements.NodesSet(e+1, self.elems[e][:])
        self.cubicWallElements.NodesSet(1, [1, 35, 2, 3, 18, 38, 19, 20])
        self.cubicWallElements.NodesSet(2, [35, 36, 3, 4, 38, 39, 20, 21])
        self.cubicWallElements.NodesSet(3, [36, 37, 4, 5, 39, 40, 21, 22])
        self.cubicWallElements.NodesSet(4, [37, 1, 5, 2, 40, 18, 22, 19])
        self.cubicWallElements.CreateFinish()
        # Set up same 16 elements mesh for pressure interpolation
        self.linWallElementMeshComponentNumber = 2
        self.linWallElements = iron.MeshElements()
        self.linWallElements.CreateStart(self.wallMesh, self.linWallElementMeshComponentNumber, self.pressureBasis)
        for e in range(4, len(self.elems)):
            self.linWallElements.NodesSet(e+1, self.elems[e][:])
        self.linWallElements.NodesSet(1, [1, 35, 2, 3, 18, 38, 19, 20])
        self.linWallElements.NodesSet(2, [35, 36, 3, 4, 38, 39, 20, 21])
        self.linWallElements.NodesSet(3, [36, 37, 4, 5, 39, 40, 21, 22])
        self.linWallElements.NodesSet(4, [37, 1, 5, 2, 40, 18, 22, 19])
        self.linWallElements.CreateFinish()
        self.wallMesh.CreateFinish()

    def _setup_wall_decomposition(self):
        # Wall mesh decomposition.
        self.wallDecompositionUserNumber = 1
        self.wallDecomposition = iron.Decomposition()
        self.wallDecomposition.CreateStart(self.wallDecompositionUserNumber, self.wallMesh)
        self.wallDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
        self.wallDecomposition.NumberOfDomainsSet(self.compNodes)
        self.wallDecomposition.CalculateFacesSet(True)
        self.wallDecomposition.CalculateLinesSet(True)
        self.wallDecomposition.CreateFinish()

    def _setup_wall_geometry(self):
        # Set up wall geometry.
        self.wallGeometricFieldUserNumber = 1
        self.wallGeometricField = iron.Field()
        self.wallGeometricField.CreateStart(self.wallGeometricFieldUserNumber, self.wallRegion)
        self.wallGeometricField.MeshDecompositionSet(self.wallDecomposition)
        self.wallGeometricField.VariableLabelSet(iron.FieldVariableTypes.U, "LVWall")
        self.wallGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
        self.wallGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
        self.wallGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
        self.wallGeometricField.ScalingTypeSet(self.scalingType)
        self.wallGeometricField.CreateFinish()
        self.wallGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
        # Create reference geometry for nodes which don't need versions.
        for i in range(0, len(self.nodes)):
            for k in range(0, 3):
                    for j in range(0, 8):
                        self.wallGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                         iron.FieldParameterSetTypes.VALUES, 1, j+1,
                                                                         i+1, k+1, self.nodes[i][k][0][j])
        # Initialise reference positions for the other three versions of node 1.
        for i in range(1, 4):
            for k in range(0, 3):
                for j in range(0, 8):
                        self.wallGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                         iron.FieldParameterSetTypes.VALUES, 1, j+1,
                                                                         i+34, k+1, self.nodes[0][k][i][j])
        # Initialise reference positions for the other three versions of node 18.
        for i in range(1, 4):
            for k in range(0, 3):
                for j in range(0, 8):
                    for v in range(1, 4):
                        self.wallGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                         iron.FieldParameterSetTypes.VALUES, 1, j+1,
                                                                         i+37, k+1, self.nodes[17][k][i][j])
        temp = numpy.concatenate(self.scaleFactors)
        self.wallGeometricField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.U, 1, temp)

    def _setup_fibre(self):
        self.fibreFieldUserNumber = 2
        self.fibreField = iron.Field()
        self.fibreField.CreateStart(self.fibreFieldUserNumber, self.wallRegion)
        self.fibreField.TypeSet(iron.FieldTypes.FIBRE)
        self.fibreField.MeshDecompositionSet(self.wallDecomposition)
        self.fibreField.GeometricFieldSet(self.wallGeometricField)
        self.fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")
        self.fibreField.ScalingTypeSet(self.scalingType)
        self.fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
        self.fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
        self.fibreField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
        self.fibreField.CreateFinish()
        # Fibre angles for all wall nodes.
        for n in range(0, len(self.fibre)):
            for c in range(0, 3):
                for deriv in range(0, 8):
                    rad_angle = self.fibre[n][c][deriv]*numpy.pi/180
                    self.fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, 1, deriv+1, n+1, c+1,
                                                             rad_angle)
        # Hanging apex nodes:
        for n in range(34, 37):
            for c in range(0, 3):
                for deriv in range(0, 8):
                    rad_angle = self.fibre[0][c][deriv]*numpy.pi/180
                    self.fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, 1, deriv+1, n+1, c+1,
                                                             rad_angle)
        for n in range(37, 40):
            for c in range(0, 3):
                for deriv in range(0, 8):
                    rad_angle = self.fibre[17][c][deriv]*numpy.pi/180
                    self.fibreField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES, 1, deriv+1, n+1, c+1,
                                                             rad_angle)
        temp = numpy.concatenate(self.scaleFactors)
        self.fibreField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.U, 1, temp)

    def _setup_cavity_mesh(self):
        # Set up cavity mesh.
        self.cavityMeshUserNumber = 2
        self.cavityMesh = iron.Mesh()
        self.cavityMesh.CreateStart(self.cavityMeshUserNumber, self.cavityRegion, self.numOfXi)
        self.cavityMesh.NumberOfComponentsSet(1)
        self.cavityMesh.NumberOfElementsSet(16)
        self.cavityNodes = iron.Nodes()
        self.cavityNodes.CreateStart(self.cavityRegion, 40)
        self.cavityNodes.CreateFinish()
        # Set up 16 elements mesh for cavity model.
        self.cavityElementUserNumber = 1
        self.cavityElements = iron.MeshElements()
        self.cavityElements.CreateStart(self.cavityMesh, self.cavityElementUserNumber, self.cavityBasis)
        self.cavityElements.NodesSet(16, [31, 19, 28, 18, 13, 10, 17, 14])
        self.cavityElements.NodesSet(15, [30, 31, 27, 28, 12, 13, 16, 17])
        self.cavityElements.NodesSet(14, [29, 30, 26, 27, 11, 12, 15, 16])
        self.cavityElements.NodesSet(13, [19, 29, 18, 26, 10, 11, 14, 15])
        self.cavityElements.NodesSet(12, [34, 20, 31, 19, 9, 6, 13, 10])
        self.cavityElements.NodesSet(11, [33, 34, 30, 31, 8, 9, 12, 13])
        self.cavityElements.NodesSet(10, [32, 33, 29, 30, 7, 8, 11, 12])
        self.cavityElements.NodesSet(9, [20, 32, 19, 29, 6, 7, 10, 11])
        self.cavityElements.NodesSet(8, [37, 21, 34, 20, 5, 2, 9, 6])
        self.cavityElements.NodesSet(7, [36, 37, 33, 34, 4, 5, 8, 9])
        self.cavityElements.NodesSet(6, [35, 36, 32, 33, 3, 4, 7, 8])
        self.cavityElements.NodesSet(5, [21, 35, 20, 32, 2, 3, 6, 7])
        # Apical elements.
        self.cavityElements.NodesSet(4, [40, 1, 37, 21, 25, 22, 5, 2])
        self.cavityElements.NodesSet(3, [39, 40, 36, 37, 24, 25, 4, 5])
        self.cavityElements.NodesSet(2, [38, 39, 35, 36, 23, 24, 3, 4])
        self.cavityElements.NodesSet(1, [1, 38,  21, 35, 22, 23, 2, 3])
        self.cavityElements.CreateFinish()
        self.cavityMesh.CreateFinish()

    def _setup_cavity_decomposition(self):
        # Cavity mesh decomposition.
        self.cavityDecompositionUserNumber = 2
        self.cavityDecomposition = iron.Decomposition()
        self.cavityDecomposition.CreateStart(self.cavityDecompositionUserNumber, self.cavityMesh)
        self.cavityDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
        self.cavityDecomposition.NumberOfDomainsSet(self.compNodes)
        self.cavityDecomposition.CalculateFacesSet(True)
        self.cavityDecomposition.CalculateLinesSet(True)
        self.cavityDecomposition.CreateFinish()

    def _setup_cavity_geometry(self):
        # Set up cavity geometry.
        self.cavityGeometricFieldUserNumber = 2
        self.cavityGeometricField = iron.Field()
        self.cavityGeometricField.CreateStart(self.cavityGeometricFieldUserNumber, self.cavityRegion)
        self.cavityGeometricField.MeshDecompositionSet(self.cavityDecomposition)
        self.cavityGeometricField.VariableLabelSet(iron.FieldVariableTypes.U, "LVCavity")
        self.cavityGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
        self.cavityGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, 1)
        self.cavityGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
        self.cavityGeometricField.ScalingTypeSet(self.scalingType)
        self.cavityGeometricField.CreateFinish()

    def _define_mechanics(self):
        [self.equationsSetClass, self.equationsSetType] = [iron.EquationsSetClasses.ELASTICITY,
                                                           iron.EquationsSetTypes.FINITE_ELASTICITY]
        if self.TCa > 0:
            self.equationsSetSubtype = iron.EquationsSetSubtypes.GUCCIONE_ACTIVECONTRACTION
        else:
            self.equationsSetSubtype = iron.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_GUCCIONE
        [self.problemClass, self.problemType] = [iron.ProblemClasses.ELASTICITY,
                                                                      iron.ProblemTypes.FINITE_ELASTICITY]

    def _setup_equations_set(self):
        self.equationsSetUserNumber = 1
        self.equationsSetFieldUserNumber = 3
        self.equationsSetField = iron.Field()
        self.equationsSet = iron.EquationsSet()
        self.equationsSet.CreateStart(self.equationsSetUserNumber, self.wallRegion, self.fibreField,
                                      [self.equationsSetClass, self.equationsSetType, self.equationsSetSubtype],
                                      self.equationsSetFieldUserNumber, self.equationsSetField)
        self.equationsSet.CreateFinish()

    def _setup_material(self):
        self.materialFieldUserNumber = 4
        self.materialField = iron.Field()
        self.equationsSet.MaterialsCreateStart(self.materialFieldUserNumber, self.materialField)
        self.materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
        self.materialField.ScalingTypeSet(self.scalingType)
        self.materialField.CreateFinish()

    def _setup_dependent(self):
        self.dependentFieldUserNumber = 5
        self.dependentField = iron.Field()
        self.equationsSet.DependentCreateStart(self.dependentFieldUserNumber, self.dependentField)
        self.dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
        self.dependentField.ScalingTypeSet(self.scalingType)
        # Set tricubic Hermite interpolation for geometry.
        for c in [1, 2, 3]:
            self.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, c, 1)
            self.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, c, 1)
        # Set trilinear interpolation for pressure.
        self.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 4, 2)
        self.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN, 4, 2)
        self.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U, 4,
                                                      iron.FieldInterpolationTypes.NODE_BASED)
        self.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN, 4,
                                                      iron.FieldInterpolationTypes.NODE_BASED)
        self.equationsSet.DependentCreateFinish()
        # Initialise dependent field with reference geometric field.
        for c in [1, 2, 3]:
            self.wallGeometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                             iron.FieldParameterSetTypes.VALUES, c,
                                                                             self.dependentField,
                                                                             iron.FieldVariableTypes.U,
                                                                             iron.FieldParameterSetTypes.VALUES, c)
            self.dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.DELUDELN,
                                                            iron.FieldParameterSetTypes.VALUES, c, 0.0)
        temp = numpy.concatenate(self.scaleFactors)
        self.dependentField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.U, 1, temp)
        self.dependentField.ParameterSetNodeScaleFactorsSet(iron.FieldVariableTypes.DELUDELN, 1, temp)
        for n in range(0, 40):
            self.dependentField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                         iron.FieldParameterSetTypes.VALUES, 1, 1, n+1, 4, 0.0)
        self.dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
        self.dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

    def _setup_independent(self):
        self.independentFieldUserNumber = 6
        self.independentField = iron.Field()
        self.equationsSet.IndependentCreateStart(self.independentFieldUserNumber, self.independentField)
        self.independentField.VariableLabelSet(iron.FieldVariableTypes.U, "Activation")
        self.equationsSet.IndependentCreateFinish()

    def _setup_deformed(self):
        self.deformedFieldUserNumber = 9
        self.deformedField = iron.Field()
        self.deformedField.CreateStart(self.deformedFieldUserNumber, self.wallRegion)
        self.deformedField.MeshDecompositionSet(self.wallDecomposition)
        self.deformedField.TypeSet(iron.FieldTypes.GEOMETRIC)
        self.deformedField.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
        self.deformedField.ScalingTypeSet(self.scalingType)
        self.deformedField.CreateFinish()

    def _setup_hydro(self):
        self.hydroFieldUserNumber = 10
        self.hydroField = iron.Field()
        self.hydroField.CreateStart(self.hydroFieldUserNumber, self.wallRegion)
        self.hydroField.MeshDecompositionSet(self.wallDecomposition)
        self.hydroField.TypeSet(iron.FieldTypes.GENERAL)
        self.hydroField.GeometricFieldSet(self.wallGeometricField)
        self.hydroField.VariableLabelSet(iron.FieldVariableTypes.U, "HydrostaticPressure")
        self.hydroField.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)
        self.hydroField.NumberOfVariablesSet(1)

        self.hydroField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 2)
        self.hydroField.CreateFinish()

    def _setup_equations(self):
        self.equations = iron.Equations()
        self.equationsSet.EquationsCreateStart(self.equations)
        self.equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
        self.equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
        self.equationsSet.EquationsCreateFinish()

    def _update_independent(self):
        self.independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                          iron.FieldParameterSetTypes.VALUES, 1, self.activeCurrent)
        for c in [2, 3]:
            self.independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,
                                                              iron.FieldParameterSetTypes.VALUES, c, 0.0)
        print 'Increase activation TCa to: '+str(self.activeCurrent), ' kPa'

    def _setup_problem(self):
        self.problemUserNumber = 1
        self.problem = iron.Problem()

        self.problem.CreateStart(self.problemUserNumber, [self.problemClass, self.problemType])
        self.problem.CreateFinish()
        self.problem.ControlLoopCreateStart()
        self.controlLoop = iron.ControlLoop()
        self.problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], self.controlLoop)
        self.problem.ControlLoopCreateFinish()
        self.linearSolver = iron.Solver()
        self.nonLinearSolver = iron.Solver()
        self.problem.SolversCreateStart()
        self.problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, self.nonLinearSolver)
        self.nonLinearSolver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
        self.nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS)
        self.nonLinearSolver.NewtonAbsoluteToleranceSet(1e-3)
        self.nonLinearSolver.NewtonSolutionToleranceSet(1e-3)
        self.nonLinearSolver.NewtonConvergenceTestTypeSet(iron.NewtonConvergenceTypes.PETSC_DEFAULT)
        self.nonLinearSolver.NewtonLinearSolverGet(self.linearSolver)
        self.nonLinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.CUBIC)
        self.nonLinearSolver.NewtonMaximumIterationsSet(100)
        self.linearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
        self.problem.SolversCreateFinish()
        self.solverEquations = iron.SolverEquations()
        self.solver = iron.Solver()
        self.problem.SolverEquationsCreateStart()
        self.problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, self.solver)
        self.solver.SolverEquationsGet(self.solverEquations)
        self.solverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
        self.solverEquations.EquationsSetAdd(self.equationsSet)
        self.problem.SolverEquationsCreateFinish()

    def _setup_boundary_conditions(self, warm_solve_toggle):
        # Basal displacement
        self.boundaryConditions = iron.BoundaryConditions()
        self.solverEquations.BoundaryConditionsCreateStart(self.boundaryConditions)
        d_fix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,
                 iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3]
        for n in self.epiBaseNodes:
            for c in [1, 2, 3]:
                if warm_solve_toggle == 0:
                    self.boundaryConditions.SetNode(self.dependentField, iron.FieldVariableTypes.U, 1,
                                                    iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, n, c,
                                                    iron.BoundaryConditionsTypes.FIXED,
                                                    self.baseBC[n-self.epiBaseNodes[0], c-1])
                elif warm_solve_toggle == 1:
                    self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1,
                                                    iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, n, c,
                                                    iron.BoundaryConditionsTypes.FIXED, 0.0)
                for d in d_fix:
                    self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, d, n, c,
                                                    iron.BoundaryConditionsTypes.FIXED, 0.0)
        # Endocardial pressure
        for n in self.endoNodes:
            if warm_solve_toggle == 0:
                self.boundaryConditions.SetNode(self.dependentField, iron.FieldVariableTypes.DELUDELN, 1, 1, n, 3,
                                                iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, self.endoPressure)
            elif warm_solve_toggle == 1:
                self.boundaryConditions.SetNode(self.dependentField, iron.FieldVariableTypes.DELUDELN, 1, 1, n, 3,
                                                iron.BoundaryConditionsTypes.PRESSURE, self.endoPressure)
        """
        for n in self.endoNodes:
            if warm_solve_toggle == 0:
                self.boundaryConditions.SetNode(self.dependentField, iron.FieldVariableTypes.DELUDELN, 1, 1, n, 3,
                                                iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, self.endoPressure)
            elif warm_solve_toggle == 1:
                self.boundaryConditions.SetNode(self.dependentField, iron.FieldVariableTypes.DELUDELN, 1, 1, n, 3,
                                                iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, self.endoPressure)
        """
        # Fix apical derivatives.
        d_fix = [iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                 iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                 iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2_S3,
                 iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3]
        for n in [1, 35, 36, 37, 18, 38, 39, 40]:
            for c in [1, 2, 3]:
                for d in d_fix:
                    self.boundaryConditions.SetNode(self.dependentField, iron.FieldVariableTypes.U, 1, d, n, c,
                                                    iron.BoundaryConditionsTypes.FIXED, 0.0)
        # Bind apical nodes together
        d_map = [iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
        for c in [1, 2, 3]:
            for d in d_map:
                self.boundaryConditions.ConstrainNodeDofsEqual(self.dependentField, iron.FieldVariableTypes.U, 1,
                                                               d, c, [1, 35, 36, 37])
                self.boundaryConditions.ConstrainNodeDofsEqual(self.dependentField, iron.FieldVariableTypes.U, 1,
                                                               d, c, [18, 38, 39, 40])
        self.solverEquations.BoundaryConditionsCreateFinish()

    def _update_cavity(self):
        self.cavityGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,
                                                          iron.FieldParameterSetTypes.VALUES)
        cavity_nodes = [22, 23, 24, 25, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        wall_nodes = [1, 35, 36, 37, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
        for n in range(0, len(wall_nodes)):
            for c in [1, 2, 3]:
                for d in [1, 2, 3, 4, 5, 6, 7, 8]:
                    value = self.dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                      iron.FieldParameterSetTypes.VALUES, 1, d,
                                                                      wall_nodes[n], c)
                    self.cavityGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                       iron.FieldParameterSetTypes.VALUES, 1, d,
                                                                       cavity_nodes[n], c, value)
                    sf = self.dependentField.ParameterSetNodeScaleFactorGet(iron.FieldVariableTypes.U, 1, d,
                                                                            wall_nodes[n], c)
                    self.cavityGeometricField.ParameterSetNodeScaleFactorSet(iron.FieldVariableTypes.U, 1, d,
                                                                             cavity_nodes[n], c, sf)
                self.cavityGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                   iron.FieldParameterSetTypes.VALUES, 1, 5,
                                                                   cavity_nodes[n], c, 0.0)
        # Evaluate coordinates of mid-cavity nodes.
        node_layers = [[14, 15, 16, 17], [10, 11, 12, 13], [6, 7, 8, 9], [2, 3, 4, 5]]
        central_nodes = [[18, 26, 27, 28], [19, 29, 30, 31], [20, 32, 33, 34], [21, 35, 36, 37]]
        for n in range(0, 4):
            for c in [1, 2, 3]:
                tot = 0
                for j in range(0, 4):
                    tot += self.dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES, 1, 1,
                                                                     node_layers[n][j], c)
                ave = tot/4
                for k in range(0, 4):
                    self.cavityGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                       iron.FieldParameterSetTypes.VALUES, 1, 1,
                                                                       central_nodes[n][k], c, ave)
                    self.cavityGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                       iron.FieldParameterSetTypes.VALUES, 1, 5,
                                                                       central_nodes[n][k], c, 0.0)

        self.cavityGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,
                                                           iron.FieldParameterSetTypes.VALUES)
        # Apical over-hanging nodes.
        for n in [1, 38, 39, 40]:
            for c in [1, 2, 3]:
                value = self.dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                  iron.FieldParameterSetTypes.VALUES, 1, 1,
                                                                  1, c)
                self.cavityGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                   iron.FieldParameterSetTypes.VALUES, 1, 1,
                                                                   n, c, value)
                self.cavityGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,
                                                                   iron.FieldParameterSetTypes.VALUES, 1, 5,
                                                                   n, c, 0.0)

    def _export_wall(self, key):
        print '\tExporting wall geometry ' + key + '...'
        for c in [1, 2, 3]:
            self.dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                         iron.FieldParameterSetTypes.VALUES, c,
                                                                         self.deformedField,
                                                                         iron.FieldVariableTypes.U,
                                                                         iron.FieldParameterSetTypes.VALUES, c)
        self.dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES, 4,
                                                                     self.hydroField, iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES, 1)

        #if key == "def":
        #    self.materialField.Destroy()
        #    self.dependentField.Destroy()
        export_field = iron.Fields()
        export_field.CreateRegion(self.wallRegion)
        export_field.NodesExport("LVWall_"+key, "FORTRAN")
        export_field.ElementsExport("LVWall_"+key, "FORTRAN")
        export_field.Finalise()

    def _export_cavity(self, key):
        print '\tExporting cavity geometry ' + key + '...'
        export_field = iron.Fields()
        export_field.CreateRegion(self.cavityRegion)
        export_field.NodesExport("LVCavity_"+key, "FORTRAN")
        export_field.ElementsExport("LVCavity_"+key, "FORTRAN")
        export_field.Finalise()

