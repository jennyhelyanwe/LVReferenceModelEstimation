import numpy as np
import os
from geometric_utils import *
from shutil import copy
import pickle
from sklearn import decomposition  # for PCA analysis.

class Geometry:
    def __init__(self, study_ids, study_frames, num_studies, synth_ds_num, bootstrap_itr):
        print 'Initialise Geometry'
        self.study_ids = study_ids
        self.study_frames = study_frames
        self.num_studies = num_studies
        self.bootstrap_itr = bootstrap_itr
        self.synth_ds_num = synth_ds_num
        # PCA properties.
        self.h_to_b_matrix = np.loadtxt(os.environ['MAIN']+'/h_to_b_matrix.txt', dtype='float')
        self.num_modes = 27
        self.scores = np.zeros((self.num_studies, self.num_modes))

        # Perform PCA
        print 'GEOMETRY: Generate data matrix...'
        #self._generate_data_matrix_displacement_ds_ref_only()
        self._get_data_matrix_ref()
        print 'GEOMETRY: Perform PCA...'
        self._pca()
        print 'GEOMETRY: Project each study to PCA to get initial guess for scores...'
        self._project_studies_onto_pca_ref()
        self.convert_scores()
        #print self.z_score
        #print self.scores

    def _frame_idx(self):
        self.idx = [0] * (int(self.tot) - int(self.ds) + 2)
        n = 0
        for i in range(int(self.ds), int(self.tot) + 1):
            self.idx[n] = i
            n += 1
        self.idx[n] = int(self.ed)

    def _frame_idx_synthetic(self):
        self.idx = []
        for i in range(0, (int(self.tot) - int(self.ds))+1):
            self.idx.append(i)

    def _generate_data_matrix_displacement_from_reference_synthetic(self):
        # Convert all diastolic frames to unit scale factors for conversion into Bezier.
        self.X = []
        self.data_info = []
        self.reference_models = np.zeros([self.num_studies, 960])
        # Loop through studies and evaluate deformation data.
        for i in range(0, self.num_studies):
            study_id = self.study_ids[i]
            study_frame = self.study_frames[i]
            self.ds, self.ed, self.es, self.tot = tuple(study_frame)
            self._frame_idx_synthetic()
            # Copy in template files
            os.system('cp '+ os.environ['REFERENCE_MODELS']+'/template_files/*.* '+os.environ['REFERENCE_MODELS'] + \
                      '/b' + str(self.bootstrap_itr) + '/')
            # Get reference models.
            #print 'Importing reference model for study '+study_id
            ref_path = os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/'
            os.chdir(ref_path)
            ref_file = ref_path + study_id
            #print ref_file
            unit_sf_file = ref_file + '_UnitSF'
            convert_arithmetic_to_unit_sf(ref_file, unit_sf_file)
            nodes = read_ipnode(unit_sf_file+ '.ipnode')
            temp = ipnode_to_dofs_format(nodes)
            temp = hermite_dofs_to_bezier_dofs(temp, self.h_to_b_matrix)
            self.reference_models[i, :] = np.reshape(temp, [960])
            # Get all other frames and calculate displacement from reference models.
            for j in range(self.synth_ds_num, len(self.idx)):
                file_name = os.environ['SYNTHETIC_DATA'] + '_no_offset/' + study_id + '/LVInflation_' + str(self.idx[j])
                unit_sf_file = file_name+'_UnitSF'
                convert_arithmetic_to_unit_sf(file_name, unit_sf_file)
                nodes = read_ipnode(unit_sf_file+'.ipnode')
                dofs = ipnode_to_dofs_format(nodes)
                dofs = hermite_dofs_to_bezier_dofs(dofs, self.h_to_b_matrix)
                dofs = np.reshape(dofs, [960])
                displacement = dofs - self.reference_models[i, :]
                self.X.append(displacement)
                self.data_info.append([study_id, self.idx[j]])
        self.X = np.array(self.X)
        pickle.dump(self.X, open('data_matrix.pl', 'w'))
        pickle.dump(self.data_info, open('data_info.pl', 'w'))
        pickle.dump(self.reference_models, open('reference_models.pl', 'w'))

    def _generate_data_matrix_displacement_ds_ref_only(self):
        # Convert all diastolic frames to unit scale factors for conversion into Bezier.
        self.X = []
        self.data_info = []
        self.reference_models = np.zeros([self.num_studies, 960])
        # Loop through studies and evaluate deformation data.
        for i in range(0, self.num_studies):
            study_id = self.study_ids[i]
            study_frame = self.study_frames[i]
            self.ds, self.ed, self.es, self.tot = tuple(study_frame)
            self._frame_idx_synthetic()
            # Copy in template files
            os.system('cp '+ os.environ['REFERENCE_MODELS']+'/template_files/*.* '+os.environ['REFERENCE_MODELS'] + \
                      '/b' + str(self.bootstrap_itr) + '/')
            # Get reference models.
            #print 'Importing reference model for study '+study_id
            ref_path = os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/'
            os.chdir(ref_path)
            ref_file = ref_path + study_id
            #print ref_file
            unit_sf_file = ref_file + '_UnitSF'
            convert_arithmetic_to_unit_sf(ref_file, unit_sf_file)
            nodes = read_ipnode(unit_sf_file+ '.ipnode')
            temp = ipnode_to_dofs_format(nodes)
            temp = hermite_dofs_to_bezier_dofs(temp, self.h_to_b_matrix)
            self.reference_models[i, :] = np.reshape(temp, [960])

            # For every study, get only the true reference state model.
            file_name = os.environ['EXTRAS'] + '/GenerateSyntheticData/' + study_id + '/ReferenceWall'
            unit_sf_file = file_name+'_UnitSF'
            convert_arithmetic_to_unit_sf(file_name, unit_sf_file)
            nodes = read_ipnode(unit_sf_file+'.ipnode')
            dofs = ipnode_to_dofs_format(nodes)
            dofs = hermite_dofs_to_bezier_dofs(dofs, self.h_to_b_matrix)
            dofs = np.reshape(dofs, [960])
            displacement = dofs - self.reference_models[i, :]
            self.X.append(displacement)
            self.data_info.append([study_id, 0])
        self.X = np.array(self.X)
        pickle.dump(self.X, open('data_matrix_ds_to_ref.pl', 'w'))
        pickle.dump(self.data_info, open('data_info_ds_to_ref.pl', 'w'))
        pickle.dump(self.reference_models, open('reference_models_ds_to_ref.pl', 'w'))

    def _get_data_matrix(self):
        self.X = pickle.load(open(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/data_matrix.pl', 'r'))
        self.data_info = pickle.load(open(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/data_info.pl', 'r'))
        self.reference_models = pickle.load(
            open(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/reference_models.pl', 'r'))

    def _get_data_matrix_ref(self):
        self.X = pickle.load(
            open(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/data_matrix_ds_to_ref.pl', 'r'))
        self.data_info = pickle.load(
            open(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/data_info_ds_to_ref.pl', 'r'))
        self.reference_models = pickle.load(
            open(os.environ['REFERENCE_MODELS'] + '/b' + str(self.bootstrap_itr) + '/reference_models_ds_to_ref.pl',
                 'r'))

    def _pca(self):
        pca = decomposition.PCA(n_components=self.num_modes)
        pca.fit(self.X)
        self.pca_mean = pca.mean_
        self.pca_components = pca.components_.T
        self.pca_variance = pca.explained_variance_
        self.pca_explained_variance = pca.explained_variance_ratio_
        self.sum_pca_explained_variance = sum(self.pca_explained_variance)
        progressive_sum_pca_explained_variance = [self.pca_explained_variance[0]]
        for i in range(1, len(self.pca_explained_variance)):
            progressive_sum_pca_explained_variance.append(progressive_sum_pca_explained_variance[i-1]+self.pca_explained_variance[i])
        #print progressive_sum_pca_explained_variance
        #print 'Sum of PCA explained variance for ' + str(self.num_modes) + ' modes:'
        #print self.sum_pca_explained_variance

    def _project_studies_onto_pca(self):
        for study_idx in range(0, 28):
            idx = self.data_info.index([self.study_ids[study_idx], self.synth_ds_num])
            Y = self.X[idx, :]
            subject_0 = Y - self.pca_mean
            self.scores[study_idx, :] = np.dot(subject_0, self.pca_components)

    def _project_studies_onto_pca_ref(self):
        for study_idx in range(0, 28):
            idx = self.data_info.index([self.study_ids[study_idx], 0])
            Y = self.X[idx, :]
            subject_0 = Y - self.pca_mean
            self.scores[study_idx, :] = np.dot(subject_0, self.pca_components)

    def _write_modes_to_file(self):
        np.savetxt('Modes.txt', self.pca_components)
        np.savetxt('PopulationWeights.txt', self.pca_explained_variance)
        np.savetxt('MeanBezierDOFs.txt', self.pca_mean)
        np.savetxt('ProjectionScores.txt', self.scores)

    def reconstruct_reference_model(self, scores, study_idx):
        print 'Reconstructing reference model. '
        study_name = self.study_ids[study_idx]
        temp = np.dot(scores, self.pca_components.T)
        reconst_mesh = temp + self.pca_mean + self.reference_models[study_idx]
        nodes = bezier_dofs_to_hermite_dofs(np.reshape(reconst_mesh, [320, 3]), self.h_to_b_matrix)
        nodes = dofs_to_ipnode_format(nodes)
        ref_file = os.environ['SIMULATIONS_OPENCMISS']+'_withstiffness/'+study_name+'/ReferenceWall'
        export_geometry(nodes, ref_file)
        convert_unit_to_arithmetic_sf(ref_file, ref_file)
        print 'Finished reconstructing reference model.'

    def convert_scores(self):
        self.std = np.std(self.scores)
        self.mean = np.average(self.scores)
        self.z_score = np.zeros(self.scores.shape)
        for study_idx in range(0, 28):
            self.z_score[study_idx, :] = (self.scores[study_idx, :] - self.mean)/self.std

    def convert_scores_back(self, z_score_temp):
        return z_score_temp * self.std + self.mean
