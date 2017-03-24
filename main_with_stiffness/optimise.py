import re
import nlopt
from simulation_opencmiss import *
from scipy.optimize import fmin_l_bfgs_b
from bkcolours import *
import numpy

class Optimise:
    def __init__(self, geometry, study_id, study_frame, c1_init, study_idx, synth_ds_num, bootstrap_itr):
        self.geometry = geometry
        self.study_id = study_id
        self.study_idx = study_idx
        self.study_frame = study_frame
        self.synth_ds_num = synth_ds_num
        self.bootstrap_itr = bootstrap_itr
        self.c1_init = c1_init
        self.work_dir = os.environ['SIMULATIONS_OPENCMISS'] + '_withstiffness/' + self.study_id + '/'
        if not os.path.exists(self.work_dir):
            print 'Making directory for study ' + self.study_id
            os.mkdir(self.work_dir)
        if not os.path.exists(self.work_dir + '/fitting_error/'):
            os.mkdir(self.work_dir + '/fitting_error/')
        os.chdir(self.work_dir)
        geom_template_files = os.environ['SIMULATIONS_OPENCMISS'] + '_withstiffness/geometric_template_files/'
        os.system('cp ' + geom_template_files + '*.* ' + self.work_dir)
        reference_model_init = os.environ['REFERENCE_MODELS'] + '_withstiffness/b' + str(self.bootstrap_itr) + '/' + self.study_id
        os.system('cp ' + reference_model_init + '.exnode ' + self.work_dir + '/ReferenceWall.exnode')
        os.system('cp ' + reference_model_init + '.exelem ' + self.work_dir + '/ReferenceWall.exelem')
        os.system('cp ' + reference_model_init + '.ipnode ' + self.work_dir + '/ReferenceWall.ipnode')
        os.system('cp ' + reference_model_init + '.ipelem ' + self.work_dir + '/ReferenceWall.ipelem')
        print 'OPTIMISE: Initialise simulation... '
        self.simulation = Simulation(self.study_id, self.study_frame, self.study_idx, self.synth_ds_num)
        self.scale = 5.6

    def optimise(self):
        # Optimisation set up.
        with open('../../EstimatedC1.txt', 'r') as f:
            data = f.readlines()
        #scores_init = [self.geometry.scores[self.study_idx][0]/100.0, self.geometry.scores[self.study_idx][1]/10.0,
        #               self.geometry.scores[self.study_idx][2]/10.0, c1_init]
        #scores_init = [(self.geometry.scores[self.study_idx][0]-150)/100.0, c1_init]
        #scores_init = [(self.geometry.scores[self.study_idx][0]-20)/100.0, c1_init*10]
        c1_init = self.c1_init
        scores_init = [self.geometry.z_score[self.study_idx][0]/self.scale, c1_init]
        print scores_init
        opt = nlopt.opt(nlopt.LN_COBYLA, 2)
        opt.set_min_objective(lambda x, grad: self.objective_function(x, self.scale))
        opt.set_lower_bounds([-100, 0.0])
        opt.set_upper_bounds([100, 100])
        opt.set_initial_step(-0.01)
        opt.set_ftol_abs(0.001)
        # Perform optimisation
        self.scores_opt = opt.optimize(scores_init)

        # Get optimal results.
        self.opt_val = opt.last_optimum_value()
        self.result = opt.last_optimize_result()
        self._display_results()
        # Double check optimised results.
        mse = self.objective_function(self.scores_opt, self.scale)
        if (mse - self.opt_val) < 1e-4:
            print 'Double checked objective function at optimised C1. '
        else:
            print bkcolours.FAIL + 'ERROR: Objective function at optimal C1: ' + str(
                mse) + ' doesn''t match what optimiser ' \
                       'gives: ' + str(self.opt_val) + bkcolours.ENDC
        # Copy over optimised reference model to next bootstrap iteration.
        print 'cp ReferenceWall.exnode ' + os.environ['REFERENCE_MODELS'] + '_withstiffness/b' + str(
            self.bootstrap_itr + 1) + '/' + self.study_id + '.exnode'
        os.system('cp ReferenceWall.exnode '+os.environ['REFERENCE_MODELS'] + '_withstiffness/b' + str(
            self.bootstrap_itr + 1)+'/'+self.study_id+'.exnode')
        os.system('cp ReferenceWall.exelem ' + os.environ['REFERENCE_MODELS'] + '_withstiffness/b' + str(
            self.bootstrap_itr + 1) + '/' + self.study_id + '.exelem')
        os.system('cp ReferenceWall.ipnode ' + os.environ['REFERENCE_MODELS'] + '_withstiffness/b' + str(
            self.bootstrap_itr + 1) + '/' + self.study_id + '.ipnode')
        os.system('cp ReferenceWall.ipelem ' + os.environ['REFERENCE_MODELS'] + '_withstiffness/b' + str(
            self.bootstrap_itr + 1) + '/' + self.study_id + '.ipelem')

    def _display_results(self):
        print bkcolours.OKGREEN + 'OPTIMISE: The optimised scores for reference model are ' + str(
            self.scores_opt[0]) + ' and c1 '+str(self.scores_opt[1])+ bkcolours.ENDC
        print bkcolours.OKGREEN + 'OPTIMISE: The optimised objective function value is ' + str(self.opt_val) + bkcolours.ENDC
        if self.result > 0:
            print bkcolours.OKGREEN + 'Optimiser log: Successful termination' + bkcolours.ENDC
            if self.result == 1:
                print bkcolours.OKGREEN + 'Optimiser log: stopped with generic success.' + bkcolours.ENDC
            elif self.result == 2:
                print bkcolours.OKGREEN + 'Optimiser log: stopped because stopval was reached.' + bkcolours.ENDC
            elif self.result == 3:
                print bkcolours.OKGREEN + 'Optimiser log: stopped because ftol_rel or ftol_abs was reached.' + bkcolours.ENDC
            elif self.result == 4:
                print bkcolours.OKGREEN + 'Optimiser log: stopped because xtol_rel or xtol_abs was reached.' + bkcolours.ENDC
            elif self.result == 5:
                print bkcolours.WARNING + 'Optimiser log: stopped because maxeval was reached.' + bkcolours.ENDC
            elif self.result == 6:
                print bkcolours.WARNING + 'Optimiser log: stopped because maxtime was reached.' + bkcolours.ENDC
        else:
            print bkcolours.FAIL + 'Optimiser log: Error code termination' + bkcolours.ENDC
            if self.result == -1:
                print bkcolours.FAIL + 'Optimiser log: stopped with generic failure.' + bkcolours.ENDC
            elif self.result == -2:
                print bkcolours.FAIL + 'Optimiser log: invalid arguments - e.g. lower bounds are bigger than upper bounds, ' \
                                      'or unknown algorithm.' + bkcolours.ENDC
            elif self.result == -3:
                print bkcolours.FAIL + 'Optimiser log: ran out of memory.' + bkcolours.ENDC
            elif self.result == -4:
                print bkcolours.FAIL + 'Optimiser log: halted due to roundoff errors limiting progress. ' + bkcolours.ENDC
            elif self.result == -5:
                print bkcolours.FAIL + 'Optimiser log: halted because of forced termination - user called nlopt_force_stop(' \
                                      'opt).' + bkcolours.ENDC

    def objective_function(self, params, scale):
        #print bkcolours.OKGREEN + 'OPTIMISE: Evaluating scores: ' + str(params[0]*100.0) + ' '+ str(params[1]*10) + ' '+ str(params[2]*10) + bkcolours.ENDC
        #scores_scaled = self.geometry.scores[self.study_idx]
        scores = self.geometry.z_score[self.study_idx]
        scores_temp = np.concatenate(([params[0]*scale], scores[1:27]))
        scores_converted = self.geometry.convert_scores_back(scores_temp)
        #scores_scaled[0] = params[0]*100
        #print scores_converted
        print bkcolours.OKGREEN + 'OPTIMISE: Evaluating 1st mode score: ' + str(params[0]) + bkcolours.ENDC
        print bkcolours.OKGREEN + 'OPTIMISE: with C1: ' + str(params[1]) + bkcolours.ENDC
        os.chdir(os.environ['SIMULATIONS_OPENCMISS'] + '_withstiffness/' + self.study_id + '/')
        # Reconstruct current guess for reference model.
        self.geometry.reconstruct_reference_model(scores_converted, self.study_idx)
        # Perform forward solve.
        mse = self.simulation.forward_solve(self.study_idx, params[1])
        return mse

    def projection_error(self, study_idx):
        # This function extracts the MSE of fitting for both endocardial and epicardial surfaces for all frames from DS+1 to
        # ED and sums them up.
        os.chdir(self.simulation.work_dir)
        mse = []
        for i in range(1, len(self.simulation.diastolic_pressure)):
            # Vector error value
            info = open('EpiError_'+str(i)+'.opdata').read()
            array_ = re.split(r'[\n\\]', info)
            temp = array_[3].split()
            num_epi = float(temp[len(temp)-1])
            temp = array_[7].split()
            epi_rmse = float(temp[len(temp)-1])
            info = open('EndoError_'+str(i)+'.opdata').read()
            array_ = re.split(r'[\n\\]', info)
            temp = array_[3].split()
            num_endo = float(temp[len(temp)-1])
            temp = array_[7].split()
            endo_rmse = float(temp[len(temp)-1])
            mse.append((epi_rmse**2*num_epi + endo_rmse**2*num_endo)/(num_endo+num_epi))
            print '\033[0;30;45m LOG: Current MSE for frame '+str(i)+' = '+str(mse[i-1])+'\033[0m\n'
        mse_tot = 0.0
        for i in range(1, len(mse)):
            mse_tot = mse_tot + mse[i-1]
        print 'LOG: Current total MSE for diastole = '+str(mse_tot)+'\n'
        return mse_tot

    def evaluate_hessian(self, scale):
        print 'OPTIMISE: Evaluating Hessian matrix...'
        self.scores_opt = [self.geometry.z_score[0][0]/scale, 2.85]
        n = len(self.scores_opt)
        h = 0.05
        ee = np.zeros([n, n])
        for i in range(0, n):
            ee[i, i] = h

        # First order derivatives
        A = np.zeros(n)
        B = np.zeros(n)
        g = np.zeros(n)
        for i in range(0, n):
            # Central difference proximation
            A[i] = self.objective_function(self.scores_opt+ee[:, i], scale)
            B[i] = self.objective_function(self.scores_opt-ee[:, i], scale)
            g[i] = (A[i] - B[i])/(2*ee[i,i])

        # Second order derivatives.
        C = np.zeros(n)
        D = np.zeros(n)
        E = np.zeros(n)
        F = np.zeros(n)
        H = np.zeros([n, n])
        for i in range(0, n):
            C = self.objective_function(self.scores_opt+ee[:, i]+ee[:, i], scale)
            E = self.objective_function(self.scores_opt, scale)
            F = self.objective_function(self.scores_opt-ee[:, i]-ee[:, i], scale)
            H[i,i] = (-C + 16*A[i] - 30*E + 16*B[i] - F)/(12*ee[i,i]*ee[i,i])
            for j in range(i+1, n):
                G = self.objective_function(self.scores_opt+ee[:, i]+ee[:, j], scale)
                I = self.objective_function(self.scores_opt+ee[:, i]-ee[:, j], scale)
                J = self.objective_function(self.scores_opt-ee[:, i] + ee[:, j], scale)
                K = self.objective_function(self.scores_opt-ee[:, i]-ee[:, j], scale)
                H[i, j] = (G - I - J + K)/(4*ee[i, i]*ee[j, j])
                H[j, i] = H[i, j]
        # Calculate determinant of Hessian matrix
        determinant = np.linalg.det(H)

        # Calculate scaled Hessian matrix
        scaled_H = np.zeros([n, n])
        for i in range(0, n):
            for j in range(0, n):
                scaled_H[i, j] = H[i, j]/np.sqrt(H[i, i]*H[j, j])

        # Write out results to text file.
        with open('Hessian_'+self.study_id + '.txt', 'w') as f:
            f.write('The Hessian matrix is:\n')
            for i in range(0, n):
                for j in range(0, n):
                    f.write(str(H[i, j]) + '\t')
                f.write('\n')
            f.write('The determinant of the Hessian matrix is:\n'+str(determinant) + '\n')
            f.write('The scaled Hessian matrix is:\n')
            for i in range(0, n):
                for j in range(0, n):
                    f.write(str(scaled_H[i, j]) + '\t')
                f.write('\n')

        print 'OPTIMISE: Scaled Hessian matrix: ' + str(scaled_H)
        print 'OPTIMISE: Determinant of Hessian: ' + str(determinant)

    def parameter_sweep(self, scale):
        #change_scores = [-40, -20, 0, 20, 40]
        #percent = [0.6, 0.8, 1.0, 1.2, 1.4]
        #d = [-1, -0.5, 0.0, 0.5, 1.0]
        d = [-0.8,  -0.7, -0.5, -0.3, -0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
        #geom_percent = percent
        #change_c1 = [-1, 0, 1, 2, 3]
        with open('../../EstimatedC1.txt', 'r') as f:
            data = f.readlines()
        c1 = float(data[self.study_idx].split()[-1])
        f = open('../../sweeps/sweep_'+self.study_id+'.txt', 'a')
        f.write('C1\tScore\tMSE\n')
        opt_score = self.geometry.z_score[self.study_idx][0]/scale
        for i in range(0, len(d)):
            score_perturbed = opt_score + d[i]
            for j in range(0, len(d)):
                print 'OPTIMISE: '+str(score_perturbed)
                c1_perturbed = c1 + d[j]
                print 'OPTIMISE: '+str(c1_perturbed)
                input = [score_perturbed, c1_perturbed]
                print 'OPTIMISE: Input: ' + str(input)
                mse = self.objective_function(input, scale)
                f.write(str(c1_perturbed)+'\t'+str(score_perturbed)+'\t'+str(numpy.sqrt(mse))+'\n')
                print 'C1: '+str(c1_perturbed)+' score: '+str(score_perturbed) + ' rmse: '+str(numpy.sqrt(mse)) + '\n'
        f.close()

    def parameter_sweep_investigation(self, scale):
        c1 = [1.225, 2.025, 0.925, 1.525]
        opt_score = [0.01543, -0.58456, 0.3154, -0.18456]
        if not os.path.exists('c1_ref_investigation/'):
            os.mkdir('c1_ref_investigation/')
        for i in range(0, len(c1)):
            print 'OPTIMISE: ref weight: '+str(opt_score[i]) + ' c1: '+str(c1[i])
            input = [opt_score[i], c1[i]]
            mse = self.objective_function(input, scale)
            os.chdir(self.study_id+'/')
            folder_name = 'c1_ref_investigation/c1_'+str(c1[i])+'_ref_'+str(opt_score[i])+'/'
            if not os.path.exists(folder_name):
                os.chdir('c1_ref_investigation/')
                os.mkdir('c1_'+str(c1[i])+'_ref_'+str(opt_score[i])+'/')
                os.chdir('../')
            os.system('cp LVCavity_* '+folder_name)
            os.system('cp LVInflation_* '+folder_name)
            os.system('cp Reference* '+folder_name)

    def sensitivity(self, scale):
        with open(os.environ['MAIN']+'/../EstimatedC1.txt', 'r') as f:
            data = f.readlines()
        c1 = float(data[self.study_idx].split()[-1])
        opt_score = self.geometry.z_score[self.study_idx][0]
        """
        input_opt = [opt_score/scale, c1]
        mse_opt = self.objective_function(input_opt)
        # Sensitivity in c1 direction:
        input_c1_plus = [opt_score/scale, c1 + 0.01]
        mse_c1_plus = self.objective_function(input_c1_plus)
        s_c1 = mse_c1_plus - mse_opt
        """
        s_c1 = 0.00160465446774
        print 'OPTIMISE: Current scale: ' + str(scale)
        # Sensitivity in w direction:
        input_score_plus = [opt_score/scale + 0.01, c1]
        mse_score_plus = self.objective_function(input_score_plus, scale)
        mse_opt = 0.00115985050418
        s_score = mse_score_plus - mse_opt
        print 'OPTIMISE: C1 sensitivity: '+ str(s_c1)
        print 'OPTIMISE: Score sensitivity: '+ str(s_score)
        print 'OPTIMISE: Scale for score: '+ str(s_c1/s_score)
        print 'OPTIMISE: Obj function: '+ str(np.abs(s_c1/s_score - 1.0))
        return np.abs(s_c1/s_score - 1.0)

    def optimise_scaling(self):
        opt_s = nlopt.opt(nlopt.LN_COBYLA, 1)
        opt_s.set_min_objective(lambda x, grad: self.sensitivity(x))
        opt_s.set_lower_bounds(0.1)
        opt_s.set_upper_bounds(100)
        opt_s.set_initial_step(0.1)
        opt_s.set_ftol_abs(0.01)
        scale = [2.5]
        # Perform optimisation
        scale_opt = opt_s.optimize(scale)

        print 'OPTIMISE: Optimal scaling for weight parameter is: ' + str(scale_opt)
