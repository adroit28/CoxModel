# CoxModel
Selection of best Cox Model using different methodologies and different selection criterion

Automatic Model Selection: Cox Proportional Hazard Model

Crux of Model
Cox_Model_Selection: Fits the model using CoxPHFitter and returns object with significant input covariates using any of the following 3 methodologies: 
	1.Backward Elimination
	2.Forward Selection
	3.Stepwise selection

All the above-mentioned models use following 3 criteria to eliminate or select any input covariate.
	1.P-Value
	2.Akaike Information Criterion
	3.Bayesian Information Criterion

Usage
Create an instance of class Cox_Model_Selection and input the following parameters:
	• data - Data on which model testing needs to be done
	• duration_col - Column with Information regarding time to event
	• event_col - Column with censor flag
	• step_size - Step size for gradient descent
	• show_progress - Boolean flag of whether to show intermediate iterations or not

Once the instance of class has been created it can be used to select best model using any of the above 3 methodologies for any of the criterion
	
	1.forward_selection(covariates, fixed_variables = [], threshold = 0.05, features_removed = [], criterion = ‘p’, logging = False)  
		• covariates – List of input variables on which model is to be fitted
		• fixed_variable – List of input variables that are to be fixed in the model selection, irrespective of criterion for model selection
		• threshold – upper limit for p-value
		• features_removed – Is used only when running stepwise regression, otherwise empty 
		• criterion – Criteria on the basis of which a particular value is selection: p-value, AIC, BIC
		• logging – Boolean flag of whether to print intermediate results or not

	2.backward_elimination(covariates, threshold, criterion, logging)

	3.stepwise_selection(covariates, threshold_BE, threshold_FS, logging, criterion)
		• threshold_BE - upper limit for p-value when removing variable using backward elimination
		• threshold_FS – upper limit fot p-value when selecting variable using forward selection
 

