# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 00:17:56 2018

@author: A1298
"""

import pandas, numpy
from lifelines import CoxPHFitter
from multiprocessing import Pool
import copy, inspect, math

min_IC = numpy.inf
sample_size = None

class CoxModelSelection:
    
    global logging_ 
    logging_ = False
    
    def __init__(self, data, duration_col, event_col, step_size = 0.5, show_progress = False):

        """
        data: dataframe to be fit into the model
        duration_col: column name with time to event values, must be integer
        event_col: column name with censor flag
        step_size: step by which to vary slope in gradient function
        show_progess: variable to show intermediate iteration or not
        """        
        
        self.data = data
        self.duration_col = duration_col
        self.event_col = event_col
        self.step_size = step_size
        self.show_progress = show_progress
        
    def univariate_analysis(self, covariates):

        """
        This computes the individual significance of input variable to the output variable in terms of p-value,
        z score.  
        
        References:
            None
        
        Parameters:
            covariates: list of input variable or features whose individual significance is to be tested
            to the output variable
            
        Returns:
            cox_model: a dataframe with coefficients, p value, z score, 95% confidence interval 
            for all input variable 
        """        

        if len(covariates) == 0:
            raise Exception, 'no input variables provided'
        
        # Empty DataFrame        
        cox_model = pandas.DataFrame()
        
        for feature in covariates:
            
            # Creating List of input covatiates to be tested
            col = [feature]
            col.extend([self.duration_col, self.event_col])
            
            # Calling fit function
            cph_ = fit(self.data, 
                       col, 
                       self.duration_col, 
                       self.event_col, 
                       self.show_progress, 
                       self.step_size)
            
            # Appending summary of all individual variables in dataframe
            cox_model = pandas.concat([cox_model, cph_.summary])
            
        return cox_model    
    
   
    def multivariate_analysis(self, covariates):    
        
        """
        This computes the relative significance of input variable to the output variable in terms of p-value,
        z score.  
        
        References:
            None
        
        Parameters:
            covariates: list of input variable or features whose relative significance is to be tested
            to the output variable
            
        Returns:
            cox_model: return CoxModel output object
            
        """        
        
        if len(covariates) == 0:
            raise Exception, 'no input variables provided'
        
        # Empty DataFrame        
        cox_model = pandas.DataFrame()
        
        # Creating List of input covatiates to be tested
        covariates.extend([self.duration_col, self.event_col])
        
        # Calling fit function
        cph_ = fit(self.data, 
                   covariates, 
                   self.duration_col, 
                   self.event_col, 
                   self.show_progress, 
                   self.step_size)
        
        # Final summary object
        cox_model = cph_.summary
        
        return cox_model


    def backward_elimination(self, covariates, threshold = 0.05, criterion = 'p', logging = False):

        """
        This computes the input covariates that are significant to the output variable on the basis
        of p-value using backward elemination methind
        
        References:
            https://cran.r-project.org/web/packages/SignifReg/SignifReg.pdf
        
        Parameters:
            covariates: list of input variable or features whose significance is to be tested
            to the output variable
            threshold: below this metric p-value will be considered significant, otherwise, insignificant
            default-0.05  
            criterion: criterion to be used while selection model - p-value, AIC, BIC
            
        Returns:
            cox_model: a dataframe with coefficients, p value, z score, 95% confidence interval 
            for significant input variable 
        """        

        if len(covariates) == 0:
            raise Exception, 'no input variables provided'
            
        if logging_: print '............Backward Elimination............'        
                
        variables_removed = []
        
        # Creating List of input covatiates to be tested
        final_features = copy.copy(covariates)
        final_features.extend([self.duration_col, self.event_col])

        run_simulation = True                
        while run_simulation:
            
            # Calling fit function
            cph_ = fit(self.data, 
                       final_features, 
                       self.duration_col, 
                       self.event_col, 
                       self.show_progress, 
                       self.step_size)
                
            cox_model_obj = cph_
                        
            if criterion == 'p':
                predictor_name = p_value(cox_model_obj, criterion, threshold)
            
            elif criterion in ['AIC', 'BIC']:
                
                if logging:
                    print '============================================================='
                    print cox_model_obj.summary
                    print '============================================================='
                    
                    if criterion == 'AIC': 
                        IC_ = AIC(cox_model_obj._log_likelihood, len(covariates) + 1)
                    elif criterion == 'BIC':
                        IC_ = BIC(cox_model_obj._log_likelihood, len(covariates) + 1, len(self.data.index))
                        
                    print str(criterion) + ' : ' + str(IC_)
                
                predictor_name, cph = information_criterion([] ,criterion, logging, self.data, 
                                                       copy.copy(final_features), self.duration_col, self.event_col, self.show_progress,
                                                       self.step_size)    
                            
            # If p-value is greater than threshold then remove the co-variate from the model
            if predictor_name is not None:
                final_features.remove(predictor_name)
                variables_removed.append(predictor_name)
                
            else:
                run_simulation = False
                
            # If the calling function is stepwise_selection then return from function the 
            # feature which is most significant
            if (getcalframe(inspect.currentframe()) == 'stepwise_selection'):
                    
                return variables_removed
                                
        return cox_model_obj
    

    def forward_selection(self, covariates, fixed_variables = [], threshold = 0.05, features_removed = [], criterion = 'p', logging = False):

        """
        This computes the input covariates that are significant to the output variable on the basis
        of p-value using forward selection methind
        
        References:
            https://www.stat.berkeley.edu/~blfang/STAT151A/STAT151A_lab09_demos.html
            http://www.biostat.jhsph.edu/~iruczins/teaching/jf/ch10.pdf
        
        Parameters:
            covariates: list of input variable or features whose significance is to be tested
            to the output variable
            fixed_variables: list of input covariates already fixed or added in the model
            features_removed: used only when using 'stepwise regression' - list of input variables removed
            in last step of backward elimination
            threshold: below this metric p-value will be considered significant, otherwise, insignificant
            default-0.05    
            criterion: criterion to be used while selection model - p-value
            
        Returns:
            cox_model: a dataframe with coefficients, p value, z score, 95% confidence interval 
            for significant input variable 
        """        

        if len(covariates) == 0:
            raise Exception, 'no input variables provided'
            
        if logging_: print '............Forward Selection............'    
        
        global sample_size
        sample_size = len(self.data.index)
                    
        # Initialazing temp variable
        features_added = copy.copy(fixed_variables)
        
        # Creating List of input covatiates to be tested
        features_added.extend([self.duration_col, self.event_col])
        final_features = features_added
        
        # Empty DataFrame        
        final_model = None
        
        # Pool object
        pool = Pool(processes = 3)
    
        run_simulation = True    
        while run_simulation:
            
            # Objects used to store intermediate results
            processes = []
            staging_output = []
            
            for covar in covariates:
                
                features_added = copy.copy(final_features)

                # List of covarites on which model is to be fit
                if (covar not in features_added) and (covar not in features_removed):
                    features_added.append(covar)
                else: 
                    continue
                
                args = (self.data,
                        copy.copy(features_added),
                        self.duration_col,
                        self.event_col,
                        self.show_progress,
                        self.step_size)
                
                # Calling fit function 
                process = pool.apply_async(fit, args)
                
                # Creating list of all the proecesses spawned
                processes.append(process)
                            
            # Storing output of all the processes        
            for process_ in processes:
                staging_output.append(process_.get())
                                
            if criterion in ['AIC', 'BIC']:
                feature_to_be_included, cph = information_criterion(staging_output, criterion, logging)
                
            elif criterion == 'p':
                feature_to_be_included, cph = p_value(staging_output, criterion, threshold)
                            
            # Removing feature corresponding to minimum p-value from the model provided
            # the p-value is less than the threshold
            if feature_to_be_included is not None:
                final_features.append(feature_to_be_included)
                features_added = final_features
                covariates.remove(feature_to_be_included)
                final_model_obj = cph
                
                if (logging) and (criterion in ['AIC', 'BIC']):
                    print '============================================================'
                    print cph.summary
                    print '============================================================'
                    
                    if (criterion == 'AIC'):
                        print 'AIC: ' + str(cph._AIC_)
                    elif (criterion == 'BIC'):
                        print 'BIC: ' + str(cph._BIC_)

            else:
                run_simulation = False
             
            # If the calling function is stepwise_selection then return from function the 
            # feature which is most significant   
            if (getcalframe(inspect.currentframe()) == 'stepwise_selection'):
                run_simulation = False
                
                if logging_:
                    
                    print 'variable selected using FS: ' + str(feature_to_be_included)
                    print '==========================================================='
                    
                    if final_model is not None: 
                        print final_model[['coef', 'p']]
                    else:
                        print cph.summary[['coef', 'p']]
                        
                    print '==========================================================='
  
                return feature_to_be_included 
            
            # If all the covariates are tested then stop the simulation 
            if len(covariates) == 0:
                run_simulation = False
           
        return final_model_obj

    def stepwise_selection(self, covariates, threshold_BE = 0.05, threshold_FS = 0.05, logging = False, criterion = 'p'):
        
        """
        This computes the input covariates that are significant to the output variable on the basis
        of p-value using both forward and backward selection
        
        References:
            http://people.stat.sfu.ca/~lockhart/richard/350/08_2/lectures/VariableSelection/web.pdf
            http://www.biostat.jhsph.edu/~iruczins/teaching/jf/ch10.pdf
            https://onlinecourses.science.psu.edu/stat501/node/329/
        
        Parameters:
            covariates: list of input variable or features whose significance is to be tested
            to the output variable
            threshold_BE: below this metric p-value will be considered significant, otherwise, insignificant
            default-0.05 for backward elimination
            threshold_BE: below this metric p-value will be considered significant, otherwise, insignificant
            default-0.05 for forward selection   
            logging: log intermediate steps on console
            criterion: criterion to be used while selection model - p-value
            
        Returns:
            cox_model: a dataframe with coefficients, p value, z score, 95% confidence interval 
            for significant input variable 
        """        
        
        features_selected = []
        run_simulation = True
        variable_removed_using_BE = []
        
        logging_ = logging
        
        while run_simulation:
        
            # Calling forward selection method
            variable_selected_using_FS = self.forward_selection(copy.copy(covariates), 
                                                                features_selected, 
                                                                threshold_FS, 
                                                                copy.copy(variable_removed_using_BE),
                                                                criterion,
                                                                logging)
            
            if variable_selected_using_FS is not None:
                features_selected.append(variable_selected_using_FS)
            else:
                run_simulation = False
                break
            
            if logging_ : print 'features selected using FS: ' + str(features_selected)
            
            # Calling backward selection method
            variable_removed_using_BE = self.backward_elimination(features_selected, threshold_BE, criterion, logging)
            
            if logging_ : print 'features removed using BE: ' + str(variable_removed_using_BE)
            
            # Remove all the variables than are eliminated from backward elimination from the list of 
            # variables selected using forward selection
            if len(variable_removed_using_BE) is not 0:
                for feature in variable_removed_using_BE:
                    if feature in features_selected:
                        features_selected.remove(feature)
            
            if len(features_selected) == len(covariates):
                run_simulation = False

            if logging_ : print 'features list for next iteration: ' + str(features_selected)
            
        features_selected.extend([self.duration_col, self.event_col])
        
        final_model = fit(self.data, 
                          features_selected, 
                          self.duration_col,
                          self.event_col, 
                          self.show_progress,
                          self.step_size)        

        return final_model.summary         


def fit(data, covariates, duration_col, event_col, show_progress, step_size):
    '''
    Fits the data and return an object with dataframe with coefficients, p value, z score, 95% 
    confidence interval for all input covariates 
    '''
    # Instantiating Cox Model Object
    cph = CoxPHFitter()
    
    # Call to fit method
    cph.fit(data[covariates], 
            duration_col = duration_col, 
            event_col = event_col, 
            show_progress = show_progress,
            step_size = step_size)   

    return cph 

def getcalframe(curframe):
    '''
    Return the calling frame to curframe - For example if function f1 is called by f2 then
    calling getcalframe from f1 will return f2 
    '''    
    # Getting calling frame 
    calframe = inspect.getouterframes(curframe, 3)
    
    return calframe[1][3]    
    
def information_criterion(staging_output, criterion, logging, *args):
    '''
    Returns the covariate with min AIC, BIC value  
    '''        

    feature_with_min_IC = None
    cph = None

    if (getcalframe(inspect.currentframe()) == 'backward_elimination'):  
        
        # Setting number of processes
        pool = Pool(processes = 4)
        
        # Getting list of input variables fed into cox model
        all_features = args[1]
        
        processes = []
        
        # Fetching list of input covariates
        covariates_ = [var for var in all_features if var not in [args[2], args[3]]]
        
        for covar in all_features:
            
            # Removing feature one by one and model is tested
            if covar not in [args[2], args[3]]:
                final_features = [var for var in all_features if var != covar] 
                
            else:
                continue
            
            args_ = (args[0],
                     copy.copy(final_features),
                     args[2],
                     args[3],
                     args[4],
                     args[5])
            
            # Calling fit function 
            process = pool.apply_async(fit, args_)
            
            # Creating list of all the proecesses spawned
            processes.append(process)
                        
        # Storing output of all the processes        
        for process_ in processes:
            staging_output.append(process_.get())
            
    global min_IC
        
    # Finding feature with max likelihood
    for cph_ in staging_output:
        
        # Number of featues: Input variables fed to the model and the intercept term
        number_of_features = 1 + len(cph_.summary.index.values)

        if criterion == 'AIC':
            IC = AIC(cph_._log_likelihood, number_of_features)
            
        elif criterion == 'BIC':
            IC = BIC(cph_._log_likelihood, number_of_features, sample_size)
        
        if (getcalframe(inspect.currentframe()) == 'backward_elimination'):
            
            covar_in_consideration = numpy.setdiff1d(covariates_, list(cph_.summary.index.values))[0]
            if logging: print str(criterion) +' on removing ' + str(covar_in_consideration) + ' is ' + str(IC)        
        
        elif (getcalframe(inspect.currentframe()) == 'forward_selection'):
            
            covar_in_consideration = cph_.summary.index.values[-1]
            if logging: print str(criterion) + ' on adding ' + str(covar_in_consideration) + ' is ' + str(IC)    
        
        # Finding min AIC value across all the models
        if IC < min_IC:
            feature_with_min_IC = covar_in_consideration
            min_IC = IC
            
            if criterion == 'AIC':
                cph_._AIC_ = IC
            elif criterion == 'BIC':
                cph_._BIC_ = IC
                
            cph = cph_
            
    return feature_with_min_IC, cph    

def p_value(staging_output, criterion, threshold):
    '''
    Returns the covariate with min p-value in case of forward selection and covatiate with
    max p-value in case of backward elimination provided they are less than the threshold
    '''            
    
    if (getcalframe(inspect.currentframe()) == 'forward_selection'): 

        min_p_value = -numpy.inf        
        for cph_ in staging_output:
            
            # Determining p value of most recent variable added to the model
            cox_model = cph_.summary
            covariate = cox_model.index.values[-1]
            p_value = cox_model.at[covariate, criterion]
            
            # Finding feature corresponding to minimum p-value
            if p_value < min_p_value:
                min_p_value = p_value
                feature_with_min_p_value = covariate
                cph = cph_
                
        if min_p_value < threshold:
            feature_with_min_p_value = feature_with_min_p_value
        else:
            feature_with_min_p_value = None
    
        return feature_with_min_p_value, cph   

    elif (getcalframe(inspect.currentframe()) == 'backward_elimination'):  
        
        # Finding max p-value and the correspoding covariate 
        cox_model = staging_output.summary
        max_p_value = cox_model[criterion].max()
        
        if max_p_value > threshold:
            predictor_name = cox_model.loc[cox_model[criterion] == max_p_value].index.values[0]
        else:
            predictor_name = None
            
        return predictor_name    
                
def AIC(log_likelihood, k):
    
    _AIC = (-2 * log_likelihood) + (2 * k)
    return _AIC   

def BIC(log_likelihood, k, n):

    _BIC = (-2 * log_likelihood) + (math.log(n) * k)
    return _BIC    