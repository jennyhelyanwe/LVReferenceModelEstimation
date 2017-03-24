__author__ = 'Jenny'
import os
from geometric_utils import *
from shutil import copy

class Simulation:
    def __init__(self, study_id, study_frame, study_idx, synth_ds_num):
        self.study_id = study_id
        self.study_frame = study_frame
        self.synth_ds_num = synth_ds_num
        self.work_dir = os.environ['SIMULATIONS'] + '/' + self.study_id + '/'
        if not os.path.exists(self.work_dir):
            print 'Making directory for study '+self.study_id
            os.mkdir(self.work_dir)
        os.chdir(self.work_dir)
        # Copy template files.
        mech_template_files = os.environ['SIMULATIONS']+'/mechanics_template_files/'
        os.system('cp '+mech_template_files+'*.* '+self.work_dir)
        geom_template_files = os.environ['SIMULATIONS']+'/geometric_template_files/'
        os.system('cp '+geom_template_files+'*.* '+self.work_dir)
        self.data_path = os.environ['SYNTHETIC_DATA']+'/'+self.study_id + '/'

    def setup_ref_cavity(self):
        dir_work = os.environ['SIMULATIONS'] + '/'+self.study_id
        os.chdir(dir_work)
        create_cavity_model('ReferenceWall.ipnode', 'ReferenceCavity.ipnode')
        os.system('cm LVHumanDSModel.com')

    def forward_solve(self, study_idx):
        print 'Solving forward solve for study '+self.study_id
        self.setup_ref_cavity()
        pressure = bc_pressure_get(self.study_id, self.study_frame)
        self.ds, self.ed, self.es, self.tot = tuple(self.study_frame)
        self.diastolic_pressure = [0] * (int(self.tot) - int(self.ds)  +2 - self.synth_ds_num)
        n = 0
        self.diastolic_pressure[n] = 0
        for i in range(int(self.ds)+self.synth_ds_num, int(self.tot)+1):
            self.diastolic_pressure[n] = pressure[i-1] - pressure[int(self.ds)-1]
            n += 1
        self.diastolic_pressure[n] = pressure[int(self.ed)-1] - pressure[int(self.ds)-1]
        print self.diastolic_pressure
        os.system('cp $MATERIAL/'+self.study_id+'.ipmate LV_CubicGuc.ipmate')
        for i in range(1, len(self.diastolic_pressure)):
            # Set pressure BC.
            print 'Setting pressure and displacement boundary conditions...'
            p_increm = self.diastolic_pressure[i] - self.diastolic_pressure[i-1]
            if i == 1:
                bc_pressure_set(p_increm, 'LV_CubicPreEpiBase_TEMPLATE.ipinit', 'LV_CubicPreEpiBase.ipinit')
            else:
                bc_pressure_set(p_increm, 'CurrentInflated_'+str(i-1)+'.ipinit', 'LV_CubicPreEpiBase.ipinit')
            # Set displacement BC
            node_idx = [31, 32, 33, 34]
            displacements = [[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0]]
            bc_displacement_set(node_idx, displacements, 'LV_CubicPreEpiBase.ipinit')

            # Copy data files to generic name
            copy(self.data_path+'/Surface_Points_Epi_'+str(i)+'.ipdata', 'Surface_Points_Epi.ipdata')
            copy(self.data_path+'/Surface_Points_Endo_'+str(i)+'.ipdata', 'Surface_Points_Endo.ipdata')

            # Solve inflation
            print 'Solving inflation...'
            os.system('cm SolveInitialInflation.com')

            # Save solutions
            print 'Saving solutions...'
            os.system('cp output/LVInflation.ipinit CurrentInflated_' + str(i) + '.ipinit')
            os.system('cp output/LVInflation.exnode LVInflation_'+str(i)+'.exnode')
            os.system('cp output/LVInflation.exelem LVInflation_'+str(i)+'.exelem')
            os.system('cp output/LVInflation.ipnode LVInflation_' + str(i) + '.ipnode')
            os.system('cp output/LVInflation.ipelem LVInflation_' + str(i) + '.ipelem')

            # Save error files
            print 'Save error files...'
            os.system('cp output_errors/EndoProjection.exdata EndoError_'+str(i)+'.exdata')
            os.system('cp output_errors/EndoProjection.opdata EndoError_'+str(i)+'.opdata')
            os.system('cp output_errors/EpiProjection.exdata EpiError_'+str(i)+'.exdata')
            os.system('cp output_errors/EpiProjection.opdata EpiError_'+str(i)+'.opdata')
        os.chdir('../')
