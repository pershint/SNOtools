#This library contains functions relevant to producing the correlation matrix
#Based on the phi coefficient.  The phi coefficient gives you a measure
#of how correlated two binary values are.
import maskbuilder as mb
import resultgetter as rg
import leakestimator as le

def getTitles(result_dict):
    '''
    Takes in the results from a bifurcation analysis run and returns what
    the title for cut1 and cut2 should be.  The order for what cut1 and
    cut2 is is dependent on how the analysis is done in ../BifurAnalysis,
    so if you change the analysis/result output logic, CHANGE IT HERE TOO!
    '''
    print(result_dict)
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


def getPhiRows(cut1_list, cut2_list, PC_list, isSymmetric):
    '''
    Takes in all the phi coefficients from each bifurcation analysis and
    returns what titles go on the x-axis, the y-axis, and the row structure
    of the phi coefficients that will go with them.
    '''
    column_titles = list(set(cut1_list))
    column_titles.sort()
    row_titles = list(set(cut2_list))
    row_titles.sort()
    if isSymmetric:
        #Column and Row titles will be the same, since if you have one set
        #of column/row pairs, you have the mirrored row/column pair value
        column_titles = list(set(column_titles + row_titles))
        column_titles.sort()
        row_titles = list(set(column_titles + row_titles))
        row_titles.sort()
    PC_rows = []
    for y in row_titles:
        row=[]
        for x in column_titles:
            cut1_indices = [i for i,j in enumerate(cut1_list) if j==x]
            cut2_indices = [i for i,j in enumerate(cut2_list) if j==y]
            #Now, get the common index; that is, the index that matches
            #This title pair in the PC_list
            #FIXME: Need to add in the case if the titles are the same, append 1 
            try:
                PC_index = list(set(cut1_indices) & set(cut2_indices))[0]
                row.append(PC_list[PC_index])
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
                        PC_index = list(set(cut1_indices) & set(cut2_indices))[0]
                        row.append(PC_list[PC_index])
                    except IndexError:
                        #No match in mirrored set. set this value as -2
                        row.append(-2)
                else:
                    #No match in mirrored set. set this value as -2
                    row.append(-2)
        PC_rows.append(row)
    return column_titles, row_titles, PC_rows

def buildTitlesAndPhiMatrix(allresult_filenames, acceptance_rates,isSymmetric):
    #The following takes bifurcation analysis results for single cuts and
    #plots out the pearson coefficients in a grid fashion
    cut1_list = []
    cut2_list = []
    PC_list = []
    for fname in allresult_filenames:
        result_dict = rg.GetResultDict(fname) #Has the file's results
        BA = le.BifurAnalysisRun(result_dict, acceptance_rates)
        PC_list.append(BA.pearson_coeff())
        cut1_title, cut2_title = getTitles(result_dict)
        cut1_list.append(cut1_title)
        cut2_list.append(cut2_title)
    #Now, if you have all the proper files to fill out a square grid, the
    #increments you should split the lists above into are the sqrt the
    #list length
    numcolumns = np.sqrt(float(len(allresult_filenames)))
    if numcolumns / int(numcolumns) != 1:
        print("It looks like you don't have all your results needed to make" + \
                "a nice nxn grid.  Just warning you...")
    
    isSymmetric = True
    column_titles, row_titles, PC_rows = getPhiRows(cut1_list, cut2_list, PC_list,isSymmetric)
    return column_titles, row_titles, PC_rows
