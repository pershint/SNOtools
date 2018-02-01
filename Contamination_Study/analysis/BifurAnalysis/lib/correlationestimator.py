#This library contains functions relevant to producing the correlation matrix
#Based on the phi coefficient.  The phi coefficient gives you a measure
#of how correlated two binary values are.
import maskbuilder as mb
import resultgetter as rg
import leakestimator as le
import playDarts as pd

from collections import OrderedDict
import numpy as np
import scipy
import scipy.linalg
import numpy.linalg

#FIXME: Initial probabilities shouldn't be shot with covariance matrix; should
#Be the lower level variables themselves that will make the cov. matrix elements
class CovarianceMatrix(object):
    def __init__(self, cov_matrix, variables):
        '''
        Class takes in a covariance matrix and probabilities of getting
        a "FAIL" for each cut in the bifurcation analysis.  Tools included
        for Cholesky decomposition or SVD decomposition, then random
        shooting pass/fails with correlations included.
        '''

        self.cov_matrix = cov_matrix
        self.variables = variables
        self.correlation_shooter_matrix = None


    def choleskydecompose(self):
        dimension = len(np.diagonal(self.cov_matrix))
        ain = scipy.array(self.cov_matrix)
        #icorr = (0.4 * np.identity(dimension))
        #print(icorr)
        #ain = ain + icorr
        eigenvals = scipy.linalg.eigvalsh(ain)
        for val in eigenvals:
            if val < 0:
                print("matrix must be positive definite to cholesky" +\
                        "decompose.  returning none")
                return none
        c = scipy.linalg.cholesky(ain, lower=True)
        #u = scipy.linalg.cholesky(ain, lower=false)
        self.correlation_shooter_matrix = c

    def svddecompose(self):
        dimension = len(np.diagonal(self.cov_matrix))
        ain = scipy.array(self.cov_matrix)
        U, V, S = numpy.linalg.svd(ain)
        self.correlation_shooter_matrix = U   #TODO: Do we only need U to random shoot probabilities?

    def shoot_corrcuts(self, numshots):
        fired_variables = []
        for shot in np.arange(numshots):
            variables_thisshot = []
            #First, shoot random numbers from a normal distribution
            fired_norms = pd.RandShoot(0,1,len(self.variables))
            #now, multiply your cholesky decomposition to add correlations
            corr_vector = self.correlation_shooter_matrix.dot(fired_norms)
            #Shoot the probabilities by multiplying by variance, adding mu
            corr_vector = corr_vector + self.initial_probabilities
            #Finally, shoot random numbers and build the "pass_fail" vector
            fired_variables.append(variables_thisshot)
        return fired_variables

class BA_Correlations(object):
    '''
    Class takes in a cut label map telling you what each entry in a passfail
    _shot array is.  Passfail_shots is an array of passfail_cutshots returns
    from a CovarianceMatrix.shoot_cuts() return.  What cuts are in each branch
    that will define the a,b,c,d boxes of the bifurcation analysis are given
    as arrays in branch1 and branch2.
    '''
    def __init__(self, cutlabel_map, branch1, branch2):
        self.cutlabel_map = cutlabel_map
        self.branch1 = branch1
        self.branch2 = branch2
        self.len_passfailarr = None
        self.a_vals = []
        self.b_vals = []
        self.c_vals = []
        self.d_vals = []

    def clear_vals(self):
        self.len_passfailarr = None
        self.a_vals = []
        self.b_vals = []
        self.c_vals = []
        self.d_vals = []


    def fillBABoxes(self,passfail_shots):
        a = 0
        b = 0
        c = 0
        d = 0
        if self.len_passfailarr is None:
            self.len_passfailarr = len(passfail_shots)
        elif self.len_passfailarr != len(passfail_shots):
            print("WARNING: You're feeding in a number of pass-fail shots " + \
                    "inconsistent with that used to fill the a,b,c,d values" + \
                    " so far.")
        for shot in passfail_shots:
            branch1_flags, branch2_flags = [], []
            branch1_fail, branch2_fail = False, False
            for j,cut in enumerate(shot):
                if self.cutlabel_map[j] in self.branch1:
                    branch1_flags.append(shot[j])
                elif self.cutlabel_map[j] in self.branch2:
                    branch2_flags.append(shot[j])
            branch1_flags = np.array(branch1_flags)
            branch2_flags = np.array(branch2_flags)
            if np.count_nonzero(branch1_flags) > 0:
                branch1_fail = True
            if np.count_nonzero(branch2_flags) > 0:
                branch2_fail = True
            if branch1_fail and branch2_fail:
                d+=1
            elif not branch1_fail and not branch2_fail:
                a+=1
            elif not branch1_fail and branch2_fail:
                b+=1
            elif branch1_fail and not branch2_fail:
                c+=1
        self.a_vals.append(a)
        self.b_vals.append(b)
        self.c_vals.append(c)
        self.d_vals.append(d)

def getTitles(result_dict):
    '''
    Takes in the results from a bifurcation analysis run and returns what
    the title for cut1 and cut2 should be.  The order for what cut1 and
    cut2 is is dependent on how the analysis is done in ../BifurAnalysis,
    so if you change the analysis/result output logic, CHANGE IT HERE TOO!
    '''
    DCmask_dict = mb.get_dcwords()
    if result_dict["DCandFit"] == 1:
        #There's both data cleaning and fits used
        if result_dict["ITR_only"]==1:
            #Cut1 is a DC mask and Cut2 is ITR
            cut2_title = "ITR"
            cut1_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask1"]))]
        elif result_dict["B14_only"]==1:
            #Cut1 is a DC mask and Cut2 is B14
            cut2_title = "Beta14"
            cut1_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask1"]))]
        else:
            #Cut1 is all DC and Cut2 is both Fitters
            cut1_title = "All DCs"
            cut2_title = "All Fits"
    else:
        #It's only either two DC cuts or the fitters
        if result_dict["DC_mask2"]==None:
            #It's two fitters
            cut1_title = "Beta14"
            cut2_title = "ITR"
        else:
            #It's two DC cuts
            cut1_title = cut1_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask1"]))]
            cut2_title = cut2_title = DCmask_dict[mb.int_to_binary_bit(int(result_dict["DC_mask2"]))]
    return cut1_title, cut2_title


def getPhiCorrRows(cut1_list, cut2_list, cut1_vars, cut2_vars, cutvar_dict, PC_list,
        isSymmetric):
    '''
    Takes in all the phi coefficients and correlation matrix elements
    from each bifurcation analysis and
    returns what titles go on the x-axis, the y-axis, and the row structure
    of the phi coefficients and correlated_list that will go with them.
    '''
    column_titles = list(set(cut1_list))
    column_titles.sort()
    #column_titles.reverse()
    row_titles = list(set(cut2_list))
    row_titles.sort()
    if isSymmetric:
        #Column and Row titles will be the same, since if you have one set
        #of column/row pairs, you have the mirrored row/column pair value
        column_titles = list(set(column_titles + row_titles))
        column_titles.sort()
        #column_titles.reverse()
        row_titles = list(set(column_titles + row_titles))
        row_titles.sort()
    PC_rows = []
    corr_rows = []
    for y in row_titles:
        PC_row=[]
        corr_row=[]
        for x in column_titles:
            if y==x:
                PC_row.append(1)
                corr_row.append(cutvar_dict[x])
                continue
            cut1_indices = [i for i,j in enumerate(cut1_list) if j==x]
            cut2_indices = [i for i,j in enumerate(cut2_list) if j==y]
            #Now, get the common index; that is, the index that matches
            #This title pair in the PC_list
            #FIXME: Need to add in the case if the titles are the same, append 1 
            try:
                index = list(set(cut1_indices) & set(cut2_indices))[0]
                PC_row.append(PC_list[index])
                corr_row.append(PC_list[index]*np.sqrt(cut1_vars[index]*cut2_vars[index]))
            except IndexError:
                #No match in data set.
                #Phi coefficient is symmetric in row/column. Try to get the
                #value from swapping cut1 and cut2
                if isSymmetric:
                    cut1_indices = [i for i,j in enumerate(cut1_list) if j==y]
                    cut2_indices = [i for i,j in enumerate(cut2_list) if j==x]
                    #Now, get the common index; that is, the index that matches
                    #This title pair in the PC_list
                    try:
                        index = list(set(cut1_indices) & set(cut2_indices))[0]
                        PC_row.append(PC_list[index])
                        corr_row.append(PC_list[index]*np.sqrt(cut1_vars[index]* \
                                cut2_vars[index]))
                    except IndexError:
                        #No match in mirrored set. set this value as -2
                        PC_row.append(-2)
                        corr_row.append(-2)
                else:
                    #No match in mirrored set. set this value as -2
                    PC_row.append(-2)
                    corr_row.append(-2)
        PC_rows.append(PC_row)
        corr_rows.append(corr_row)
    return column_titles, row_titles, PC_rows, corr_rows


def buildTitlesAndPhiCorrMatrices(allresult_filenames, acceptance_rates,isSymmetric):
    #The following takes bifurcation analysis results for single cuts and
    #plots out the pearson coefficients in a grid fashion
    cut1_list = []
    cut2_list = []
    cut1_variances = []
    cut2_variances = []
    cutvar_dict = {}
    PC_list = []
    correl_list = []
    #Get the bifurcation analysis results for each file and
    #Get the information relevant for the correlation and covariance matrix
    for fname in allresult_filenames:
        result_dict = rg.GetResultDict(fname) #Has the file's results
        BA = le.BifurAnalysisRun(result_dict, acceptance_rates)
        PC_list.append(BA.pearson_coeff())
        cut1_variances.append(BA.V_cut1)
        cut2_variances.append(BA.V_cut2)
        cut1_title, cut2_title = getTitles(result_dict)
        cut1_list.append(cut1_title)
        cut2_list.append(cut2_title)
        cutvar_dict[cut1_title] = BA.V_cut1
        cutvar_dict[cut2_title] = BA.V_cut2

    #Now, if you have all the proper files to fill out a square grid, the
    #increments you should split the lists above into are the sqrt the
    #list length
    numcolumns = np.sqrt(float(len(allresult_filenames)))
    if numcolumns / int(numcolumns) != 1:
        print("It looks like you don't have all your results needed to make" + \
                "a nice nxn grid.  Just warning you...")
    column_titles, row_titles, PC_rows, cov_rows = \
            getPhiCorrRows(cut1_list, cut2_list, cut1_variances, 
            cut2_variances, cutvar_dict, PC_list,isSymmetric)
    print("COVARIANCE MATRIX: " + str(cov_rows))
    return column_titles, row_titles, PC_rows, cov_rows
