import re
import nlopt
from simulation_opencmiss import *
from scipy.optimize import fmin_l_bfgs_b
from bkcolours import *


class Optimise:
    def __init__(self, geometry, study_id, study_frame, study_idx, synth_ds_num, bootstrap_itr):
        self.geometry = geometry
        self.study_id = study_id
        self.study_idx = study_idx
        self.study_frame = study_frame
        self.synth_ds_num = synth_ds_num
        self.bootstrap_itr = bootstrap_itr
        self.work_dir = os.environ['SIMULATIONS_OPENCMISS'] + '/' + self.study_id + '/'
        if not os.path.exists(self.work_dir):
            print 'Making directory for study ' + self.study_id
            os.mkdir(self.work_dir)
        if not os.path.exists(self.work_dir + '/fitting_error/'):
            os.mkdir(self.work_dir + '/fitting_error/')
        os.chdir(self.work_dir)
        geom_template_files = os.environ['SIMULATIONS_OPENCMISS'] + '/geometric_template_files/'
        os.system('cp ' + geom_template_files + '*.* ' + self.work_dir)
        reference_model_init = os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/' + self.study_id
        os.system('cp ' + reference_model_init + '.exnode ' + self.work_dir + '/ReferenceWall.exnode')
        os.system('cp ' + reference_model_init + '.exelem ' + self.work_dir + '/ReferenceWall.exelem')
        os.system('cp ' + reference_model_init + '.ipnode ' + self.work_dir + '/ReferenceWall.ipnode')
        os.system('cp ' + reference_model_init + '.ipelem ' + self.work_dir + '/ReferenceWall.ipelem')
        print 'OPTIMISE: Initialise simulation... '
        self.simulation = Simulation(self.study_id, self.study_frame, self.study_idx, self.synth_ds_num)

    def optimise(self):
        # Optimisation set up.
        scores_init = [(self.geometry.scores[self.study_idx][0]-30)/100.0]
        print scores_init
        #scores_init = [self.geometry.scores[self.study_idx][0]/100.0, self.geometry.scores[self.study_idx][1]/10.0,
        #               self.geometry.scores[self.study_idx][2]/10.0]
        """
        opt = nlopt.opt(nlopt.LN_COBYLA, 1)
        opt.set_min_objective(lambda x, grad: self.objective_function(x))
        opt.set_lower_bounds(-100)
        opt.set_upper_bounds(100)
        opt.set_initial_step(-0.1)
        opt.set_ftol_abs(0.1)
        opt.set_ftol_rel(0.1)
        opt.set_maxeval(300)
        """
        opt = nlopt.opt(nlopt.LD_SLSQP, 1)
        opt.set_min_objective(lambda x, grad: self.objective_function(x))
        opt.set_lower_bounds(-100)
        opt.set_upper_bounds(100)
        opt.set_initial_step(0.1)
        opt.set_ftol_abs(1e-5)
        opt.set_ftol_rel(1e-5)
        opt.set_maxeval(300)
        # Perform optimisation
        self.scores_opt = opt.optimize(scores_init)

        # Get optimal results.
        self.opt_val = opt.last_optimum_value()
        self.result = opt.last_optimize_result()
        self._display_results()
        # Double check optimised results.
        mse = self.objective_function(self.scores_opt)
        if (mse - self.opt_val) < 1e-4:
            print 'Double checked objective function at optimised C1. '
        else:
            print bkcolours.FAIL + 'ERROR: Objective function at optimal C1: ' + str(
                mse) + ' doesn''t match what optimiser ' \
                       'gives: ' + str(self.opt_val) + bkcolours.ENDC
        # Copy over optimised reference model to next bootstrap iteration.
        os.system('cp ReferenceWall.exnode '+os.environ['REFERENCE_MODELS'] + '/b' + str(
            self.bootstrap_itr+1)+'/'+self.study_id+'.exnode')
        os.system('cp ReferenceWall.exelem ' + os.environ['REFERENCE_MODELS'] + '/b' + str(
            self.bootstrap_itr + 1) + '/' + self.study_id + '.exelem')
        os.system('cp ReferenceWall.ipnode ' + os.environ['REFERENCE_MODELS'] + '/b' + str(
            self.bootstrap_itr + 1) + '/' + self.study_id + '.ipnode')
        os.system('cp ReferenceWall.ipelem ' + os.environ['REFERENCE_MODELS'] + '/b' + str(
            self.bootstrap_itr + 1) + '/' + self.study_id + '.ipelem')

    def _display_results(self):
        print bkcolours.OKGREEN + 'OPTIMISE: The optimised scores for reference model are ' + str(
            self.scores_opt[0] * 100)  + bkcolours.ENDC
        print bkcolours.OKGREEN + 'OPTIMISE: The optimised objective function value is ' + str(self.opt_val) + bkcolours.ENDC
        if self.result > 0:
            print bkcolours.OKGREEN + 'OPTIMISE: Successful termination' + bkcolours.ENDC
            if self.result == 1:
                print bkcolours.OKGREEN + 'OPTIMISE: stopped with generic success.' + bkcolours.ENDC
            elif self.result == 2:
                print bkcolours.OKGREEN + 'OPTIMISE: stopped because stopval was reached.' + bkcolours.ENDC
            elif self.result == 3:
                print bkcolours.OKGREEN + 'OPTIMISE: stopped because ftol_rel or ftol_abs was reached.' + bkcolours.ENDC
            elif self.result == 4:
                print bkcolours.OKGREEN + 'OPTIMISE: stopped because xtol_rel or xtol_abs was reached.' + bkcolours.ENDC
            elif self.result == 5:
                print bkcolours.WARNING + 'OPTIMISE: stopped because maxeval was reached.' + bkcolours.ENDC
            elif self.result == 6:
                print bkcolours.WARNING + 'OPTIMISE: stopped because maxtime was reached.' + bkcolours.ENDC
        else:
            print bkcolours.FAIL + 'OPTIMISE: Error code termination' + bkcolours.ENDC
            if self.result == -1:
                print bkcolours.FAIL + 'OPTIMISE: stopped with generic failure.' + bkcolours.ENDC
            elif self.result == -2:
                print bkcolours.FAIL + 'OPTIMISE: invalid arguments - e.g. lower bounds are bigger than upper bounds, ' \
                                      'or unknown algorithm.' + bkcolours.ENDC
            elif self.result == -3:
                print bkcolours.FAIL + 'OPTIMISE: ran out of memory.' + bkcolours.ENDC
            elif self.result == -4:
                print bkcolours.FAIL + 'OPTIMISE: halted due to roundoff errors limiting progress. ' + bkcolours.ENDC
            elif self.result == -5:
                print bkcolours.FAIL + 'OPTIMISE: halted because of forced termination - user called nlopt_force_stop(' \
                                      'opt).' + bkcolours.ENDC

    def objective_function(self, score):
        print bkcolours.OKGREEN + 'OPTIMISE: Evaluating 1st mode score: ' + str(score*100.0) + bkcolours.ENDC
        os.chdir(os.environ['SIMULATIONS_OPENCMISS'] + '/' + self.study_id + '/')
        scores_scaled = self.geometry.scores[self.study_idx]
        scores_scaled[0] = score*100
        print scores_scaled

        # Reconstruct current guess for reference model.
        self.geometry.reconstruct_reference_model(scores_scaled, self.study_idx)

        # Perform forward solve.
        mse = self.simulation.forward_solve(self.study_idx)
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

    def evaluate_hessian(self):
        print 'OPTIMISE: Evaluating Hessian matrix...'
        n = len(self.scores_opt)
        h = 0.1
        ee = np.zeros([n, n])
        for i in range(0, n):
            ee[i, i] = h

        # First order derivatives
        A = np.zeros(n)
        B = np.zeros(n)
        g = np.zeros(n)
        for i in range(0, n):
            # Central difference proximation
            A[i] = self.objective_function(self.scores_opt+ee[:, i])
            B[i] = self.objective_function(self.scores_opt-ee[:, i])
            g[i] = (A[i] - B[i])/(2*ee[i,i])

        # Second order derivatives.
        C = np.zeros(n)
        D = np.zeros(n)
        E = np.zeros(n)
        F = np.zeros(n)
        H = np.zeros([n, n])
        for i in range(0, n):
            C = self.objective_function(self.scores_opt+ee[:, i]+ee[:, i])
            E = self.objective_function(self.scores_opt)
            F = self.objective_function(self.scores_opt-ee[:, i]-ee[:, i])
            H[i,i] = (-C + 16*A[i] - 30*E + 16*B[i] - F)/(12*ee[i,i]*ee[i,i])
            for j in range(i+1, n):
                G = self.objective_function(self.scores_opt+ee[:, i]+ee[:, j])
                I = self.objective_function(self.scores_opt+ee[:, i]-ee[:, j])
                J = self.objective_function(self.scores_opt-ee[:, i] + ee[:, j])
                K = self.objective_function(self.scores_opt-ee[:, i]-ee[:, j])
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






