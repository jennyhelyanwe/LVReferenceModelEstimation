__author__ = 'Jenny'
import os
from geometric_utils import *
from shutil import copy
import LVFE
from bkcolours import *

class Simulation:
    def __init__(self, study_id, study_frame, study_idx, synth_ds_num):
        self.study_id = study_id
        self.study_frame = study_frame
        self.synth_ds_num = synth_ds_num
        self.data_path = os.environ['SYNTHETIC_DATA']+'_no_offset/'+self.study_id + '/'
        self.LV = LVFE.LV()


    def forward_solve(self, study_idx):
        with open('../../EstimatedC1.txt', 'r') as f:
            data = f.readlines()
        c1 = data[study_idx].split()[-1]
        # Set up OPENCMISS model.
        print 'Set up LV model'
        self.LV.setup()
        print 'LV.def_material'
        self.LV.def_material(float(c1))
        print 'Solving forward solve for study '+self.study_id
        pressure = bc_pressure_get(self.study_id, self.study_frame)
        self.ds, self.ed, self.es, self.tot = tuple(self.study_frame)
        self.diastolic_pressure = [0] * (int(self.tot) - int(self.ds) + 2 - self.synth_ds_num)
        n = 0
        self.diastolic_pressure[n] = 0
        synth_frames = [self.synth_ds_num]
        for i in range(int(self.ds)+self.synth_ds_num, int(self.tot)+1):
            self.diastolic_pressure[n] = pressure[i-1]
            n += 1
            synth_frames.append(self.synth_ds_num + n)
        self.diastolic_pressure[n] = pressure[int(self.ed)-1]
        mse = []
        print self.diastolic_pressure
        for i in range(0, len(self.diastolic_pressure)):
            # Set pressure BC.
            print 'Setting pressure and displacement boundary conditions...'
            if i == 0:
                self.LV.set_pressure(self.diastolic_pressure[i], 0.0)
            else:
                self.LV.set_pressure(self.diastolic_pressure[i], self.diastolic_pressure[i-1])
            self.LV.set_base_displacement_zero()
            if i > 0:
                self.LV.initialise_dependent(i-1)
            self.LV.solve(0)

            # Copy forward solve models from generic names.
            copy('LVWall_def.part0.exnode', 'LVInflation_' + str(i) + '.exnode')
            copy('LVWall_def.part0.exelem', 'LVInflation_' + str(i) + '.exelem')
            copy('LVCavity_def.part0.exnode', 'LVCavity_' + str(i) + '.exnode')
            copy('LVCavity_def.part0.exelem', 'LVCavity_' + str(i) + '.exelem')

            # Copy data files to generic name
            copy(self.data_path+'Surface_Points_Epi_'+str(synth_frames[i])+'.ipdata', 'Surface_Points_Epi.ipdata')
            copy(self.data_path+'Surface_Points_Endo_'+str(synth_frames[i])+'.ipdata', 'Surface_Points_Endo.ipdata')
            print 'Surface_Points_Epi_'+str(synth_frames[i])+'.ipdata'
            mse_current_frame = self.LV.projection('Surface_Points_Epi.ipdata', 'Surface_Points_Endo.ipdata', i)
            mse.append(mse_current_frame)

        mse_tot = 0
        for i in range(0, len(mse)):
            mse_tot = mse_tot + mse[i]
        os.chdir('../')
        print bkcolours.OKBLUE + 'OPTIMISE: Total MSE for current guess is ' + str(mse_tot) + bkcolours.ENDC
        #self.LV.finalise()
        self.LV.finalise()
        return mse_tot
