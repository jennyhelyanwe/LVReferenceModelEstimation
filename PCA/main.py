
import numpy as np
import os
import pickle
from utils import *
from sklearn import decomposition  # for PCA analysis.


class PCA:
    def __init__(self):
        self.study_ids = []
        self.study_frames = []
        self.num_studies = 0
        self.h_to_b_matrix = np.loadtxt('h_to_b_matrix.txt', dtype='float')
        # PCA properties.
        self.num_modes = 10

    def get_study_ids(self):
        # Get list of studies
        filename = os.environ['PARAM_ESTIMATION'] + '/NYStFranFrameNumber_UsedForAnalysis.txt'
        f = open(filename, 'r')
        study_info = f.readline()
        num_studies = 0
        while len(study_info) != 0:  # Reaching the end of the file
            self.study_ids.append(study_info.split()[0])
            self.study_frames.append(study_info.split()[1:5])
            study_info = f.readline()
            self.num_studies = len(self.study_ids)

    def generate_data_matrix(self):
        # Adapted from code by Mahyar Osanlouy 1 Nov 2016
        if not os.path.exists('Models_UnitSF'):
            os.mkdir('Models_UnitSF')
        # Get first study and use as reference for procrustean transformation
        study_id = self.study_ids[0]
        study_frame = self.study_frames[0]
        ds, ed, es, tot = tuple(study_frame)
        filename_no_extension = os.environ['GEOM_DATA'] + study_id + '/Passive/' + study_id + '_' + str(ds)
        if not os.path.exists('Models_UnitSF/DSWall_' + study_id + '.ipnode'):
            convert_arithmetic_to_unit_sf(filename_no_extension, 'Models_UnitSF/DSWall_' + study_id)
        nodes = read_ipnode('Models_UnitSF/DSWall_' + study_id + '.ipnode')
        reference_model = ipnode_to_dofs_format(nodes)
        reference_model = hermite_dofs_to_bezier_dofs(reference_model, self.h_to_b_matrix)
        self.X = np.zeros((self.num_studies, 960))
        self.tform = []
        for i in range(0, self.num_studies):
            study_id = self.study_ids[i]
            study_frame = self.study_frames[i]
            ds, ed, es, tot = tuple(study_frame)
            filename_no_extension = os.environ['GEOM_DATA'] + study_id + '/Passive/' + study_id + '_' + str(ds)
            if not os.path.exists('Models_UnitSF/DSWall_' + study_id + '.ipnode'):
                convert_arithmetic_to_unit_sf(filename_no_extension, 'Models_UnitSF/DSWall_' + study_id)
            nodes = read_ipnode('Models_UnitSF/DSWall_' + study_id + '.ipnode')

            # Convert mesh to Bezier
            dofs = ipnode_to_dofs_format(nodes)
            dofs = hermite_dofs_to_bezier_dofs(dofs, self.h_to_b_matrix)
            d, Z, tform_0 = procrustes(reference_model, dofs)
            self.tform.append(tform_0)
            #scaled_dofs = translated_dofs * tform_0['scale'] Comment out scaling for better PCA.
            rotated_dofs = np.dot(dofs,tform_0['rotation'])
            translated_dofs = rotated_dofs + np.array(tform_0['translation'])
            self.X[i, :] = np.reshape(translated_dofs, [960])

    def pca(self):
        pca = decomposition.PCA(n_components=self.num_modes)
        pca.fit(self.X)
        self.pca_mean = pca.mean_
        self.pca_components = pca.components_.T
        self.pca_variance = pca.explained_variance_
        self.pca_explained_variance = pca.explained_variance_ratio_
        self.sum_pca_explained_variance = sum(self.pca_explained_variance)
        #print self.pca_explained_variance
        progressive_sum_pca_explained_variance = [self.pca_explained_variance[0]]
        for i in range(1, len(self.pca_explained_variance)):
            progressive_sum_pca_explained_variance.append(progressive_sum_pca_explained_variance[i-1]+self.pca_explained_variance[i])
        #print progressive_sum_pca_explained_variance
        print 'Sum of PCA explained variance for ' + str(self.num_modes) + ' modes:'
        print self.sum_pca_explained_variance

    def convert_scores(self, scores, SD, mean):
        sd = np.std(scores)
        mean = np.average(scores)
        self.z_score = []
        for i in range(len(scores)):
            self.z_score.append((scores[i] - mean[i]) / SD[i])

    def compare_reconstruction(self, study_idx):
        study_name = self.study_ids[study_idx]
        print 'Project study '+study_name+' onto pca'
        Y = self.X[0, :]
        subject_0 = Y - self.pca_mean
        score_0 = np.dot(subject_0, self.pca_components)

        # Reconstruct model
        temp = np.dot(score_0,self.pca_components.T)
        reconst_mesh = temp + self.pca_mean

        # Reverse General Procrustes Alignment
        reconst_mesh = procrustes_reverse(reconst_mesh, self.tform[0])

        nodes = bezier_dofs_to_hermite_dofs(np.reshape(reconst_mesh, [320, 3]), self.h_to_b_matrix)
        nodes = dofs_to_ipnode_format(nodes)
        export_geometry(nodes, 'ReconstructionModels/Reconstructed_STF_01')

        # Evaluate RMSE of projection
        print 'Evaluate RMSE...'
        evaluate_reconstruction_error(study_name)
        with open('ReconstructionModels/'+study_name+'_Endo_Error.opdata') as f:
            endo_error = f.readlines()
        with open('ReconstructionModels/'+study_name+'_Epi_Error.opdata') as f:
            epi_error = f.readlines()
        print 'Epicardial RMSE of reconstruction is:'
        print float(endo_error[7].split()[-1])
        print 'Endocardial RMSE reconstruction is:'
        print float(epi_error[7].split()[-1])

    def write_modes_to_file(self):
        np.savetxt('Modes.txt', self.pca_components)
        np.savetxt('PopulationWeights.txt', self.pca_explained_variance)
        np.savetxt('MeanBezierDOFs.txt', self.pca_mean)
        with open('GPA.p', 'wb') as f:
            pickle.dump(self.tform, f)


def main():
    analysis = PCA()

    # Get list of 28 subjects in study and their landmark frame numbers.
    analysis.get_study_ids()

    # Generate data matrix X by concatenating all 28 studies after General Procrustes Alignment (GPA) (excluding.
    analysis.generate_data_matrix()

    # Perform principal component analysis on data matrix X.
    analysis.pca()

    # TEST: Compare reconstructed model using 5 modes to original
    analysis.compare_reconstruction(0)

    # Write modes and weight to text file.
    #analysis.write_modes_to_file()


if __name__ == '__main__':
    main()
