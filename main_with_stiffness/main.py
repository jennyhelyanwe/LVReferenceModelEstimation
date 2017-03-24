import os, sys
from shutil import copy
from geometry import *
from optimise import *
# This script handles the overarching bootstrap method for estimating reference model using PCA dimension reduction.


class Bootstrap:
    def __init__(self, study_idx, log):
        self.study_ids = []
        self.study_frames = []
        self.c1_init = []
        self.num_studies = 0
        self.bootstrap_itr = 1
        self.synth_ds_num = 0
        # Get study IDs and study frames.
        self._get_frames()

        # Set up log file.
        if log == 1:
            so = se = open(os.environ['DIAGNOSTICS'] + '/Output_withstiffness_'+str(study_idx)+'.log', 'w', 0)
            sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
            os.dup2(so.fileno(), sys.stdout.fileno())
            os.dup2(se.fileno(), sys.stderr.fileno())
        # Set up bootstrap method by copying DS models as reference models for initial guess.
        print 'Copy over synthetic DS frames as initial guess'
        self._synthetic_initial_setup()
        self.geometry = Geometry(self.study_ids, self.study_frames, self.num_studies, self.synth_ds_num, self.bootstrap_itr)
        self.iterations = 3
        # Run first iteration of bootstrap, set up optimisation for reference model.
        #for study_idx in range(25, 28):
        print 'Optimise reference model...'
        self.optimise = Optimise(self.geometry, self.study_ids[study_idx], self.study_frames[study_idx], self.c1_init[study_idx], study_idx,
                                 self.synth_ds_num, self.bootstrap_itr)
        #self.optimise.optimise()
        #self.optimise.parameter_sweep(self.optimise.scale)
        self.optimise.parameter_sweep_investigation(self.optimise.scale)
        #self.optimise.optimise_scaling()
        #quit()
        #ref_error = self.ref_error_evaluate(study_idx)
        #ds_ref_error = self.ds_to_true_ref_error(study_idx)
        #print 'OPTIMISE: Error between estimated and true reference models is: ' + str(ref_error) + ' mm RMSE.'
        #print 'OPTIMISE: Error between synth DS model and true reference models is: ' + str(ds_ref_error) + ' mm RMSE.'
        #self.optimise.evaluate_hessian(self.optimise.scale)

    def _get_frames(self):
        # This function gets the specified study ID and important frame numbers.
        filename = os.environ['HEART_FAILURE_ROOT'] + '/NYStFranFrameNumber_UsedForAnalysis.txt'
        f = open(filename, 'r')
        study_info = f.readline()
        while len(study_info) != 0:  # Reaching the end of the file
            self.study_ids.append(study_info.split()[0])
            self.study_frames.append(study_info.split()[1:5])
            self.c1_init.append(float(study_info.split()[5]))
            study_info = f.readline()
            self.num_studies = len(self.study_ids)

    def _synthetic_initial_setup(self):
        if not os.path.exists(os.environ['REFERENCE_MODELS']+'/b'+str(self.bootstrap_itr)):
            os.mkdir(os.environ['REFERENCE_MODELS']+'/b'+str(self.bootstrap_itr))
        for i in range(0, self.num_studies):
            # Use first frame after reference model inflation as synthetic diastasis model.
            # In non-synthetic analysis, I would use the actual DS model geometries from CIM.
            from_file = os.environ['SYNTHETIC_DATA'] + '_no_offset/'+self.study_ids[i] + '/LVInflation_'+str(self.synth_ds_num)+'.ipnode'
            copy(from_file, os.environ['REFERENCE_MODELS']+'/b'+str(self.bootstrap_itr)+'/'+self.study_ids[i]+'.ipnode')
            from_file = os.environ['SYNTHETIC_DATA'] + '_no_offset/'+self.study_ids[i] + '/LVInflation_'+str(self.synth_ds_num)+'.ipelem'
            copy(from_file,
                 os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/' + self.study_ids[i] + '.ipelem')
            from_file = os.environ['SYNTHETIC_DATA'] + '_no_offset/' + self.study_ids[i] + '/LVInflation_'+str(self.synth_ds_num)+'.exnode'
            copy(from_file,
                 os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/' + self.study_ids[i] + '.exnode')
            from_file = os.environ['SYNTHETIC_DATA'] + '_no_offset/' + self.study_ids[i] + '/LVInflation_'+str(self.synth_ds_num)+'.exelem'
            copy(from_file,
                 os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/' + self.study_ids[i] + '.exelem')
        if not os.path.exists(os.environ['REFERENCE_MODELS']+'/b'+str(self.bootstrap_itr+1)):
            os.mkdir(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr+1))

    def ref_error_evaluate(self, study_idx):
        # Create surface points from estimated reference model, and project onto true reference models to calculate error.
        file = os.environ['SIMULATIONS_OPENCMISS']+'_withstiffness/'+self.study_ids[study_idx]+'/ReferenceWall'
        true_ref = os.environ['EXTRAS'] + '/GenerateSyntheticData/'+self.study_ids[study_idx]+'/ReferenceWall'
        os.chdir(os.environ['SIMULATIONS_OPENCMISS']+'_withstiffness/'+self.study_ids[study_idx])
        with open('RefErrorCalc.com', 'r') as f:
            data = f.readlines()

        data[5] = 'fem def node;r;'+file+'\n'
        data[6] = 'fem def elem;r;'+file+'\n'
        data[28] = 'fem def node;r;'+true_ref+'\n'
        data[29] = 'fem def elem;r;'+true_ref+'\n'
        with open('RefErrorCalc.com', 'w') as f:
            f.writelines(data)

        os.system('cm RefErrorCalc.com')

        copy('RefSurfaceErrorEndo.opdata', self.study_ids[study_idx]+'_RefSurfaceErrorEndo.opdata')
        copy('RefSurfaceErrorEndo.exdata', self.study_ids[study_idx] + '_RefSurfaceErrorEndo.exdata')
        copy('RefSurfaceErrorEpi.opdata', self.study_ids[study_idx] + '_RefSurfaceErrorEpi.opdata')
        copy('RefSurfaceErrorEpi.exdata', self.study_ids[study_idx] + '_RefSurfaceErrorEpi.exdata')

        with open('RefSurfaceErrorEndo.opdata', 'r') as f:
            data = f.readlines()
        endo_error = float(data[7].split()[-1])
        num_endo_data = float(data[3].split()[-1])
        with open('RefSurfaceErrorEpi.opdata', 'r') as f:
            data = f.readlines()
        epi_error = float(data[7].split()[-1])
        num_epi_data = float(data[3].split()[-1])
        ref_error = np.sqrt((endo_error**2*num_endo_data + epi_error**2*num_epi_data)/(num_endo_data + num_epi_data))
        return ref_error

    def ds_to_true_ref_error(self, study_idx):
        # Create surface points from estimated reference model, and project onto true reference models to calculate error.
        file = os.environ['SYNTHETIC_DATA'] +'_no_offset/'+self.study_ids[study_idx] + '/LVInflation_' + str(self.synth_ds_num)
        true_ref = os.environ['EXTRAS'] + '/GenerateSyntheticData/' + self.study_ids[study_idx] + '/ReferenceWall'
        os.chdir(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr))
        with open('RefErrorCalc.com', 'r') as f:
            data = f.readlines()

        data[5] = 'fem def node;r;' + file + '\n'
        data[6] = 'fem def elem;r;' + file + '\n'
        data[28] = 'fem def node;r;' + true_ref + '\n'
        data[29] = 'fem def elem;r;' + true_ref + '\n'
        with open('RefErrorCalc.com', 'w') as f:
            f.writelines(data)

        os.system('cm RefErrorCalc.com')

        copy('RefSurfaceErrorEndo.opdata', self.study_ids[study_idx] + '_RefSurfaceErrorEndo.opdata')
        copy('RefSurfaceErrorEndo.exdata', self.study_ids[study_idx] + '_RefSurfaceErrorEndo.exdata')
        copy('RefSurfaceErrorEpi.opdata', self.study_ids[study_idx] + '_RefSurfaceErrorEpi.opdata')
        copy('RefSurfaceErrorEpi.exdata', self.study_ids[study_idx] + '_RefSurfaceErrorEpi.exdata')

        with open('RefSurfaceErrorEndo.opdata', 'r') as f:
            data = f.readlines()
        endo_error = float(data[7].split()[-1])
        num_endo_data = float(data[3].split()[-1])
        with open('RefSurfaceErrorEpi.opdata', 'r') as f:
            data = f.readlines()
        epi_error = float(data[7].split()[-1])
        num_epi_data = float(data[3].split()[-1])
        ref_error = np.sqrt(
            (endo_error ** 2 * num_endo_data + epi_error ** 2 * num_epi_data) / (num_endo_data + num_epi_data))
        return ref_error

study_idx = int(sys.argv[1])
log = int(sys.argv[2])
analysis = Bootstrap(study_idx, log)
