"""
The analysis of the EOF and MCA

.. seealso::
    `xeofs`: https://xeofs.readthedocs.io/en/latest/
"""
import xarray as xr
import xeofs
import pickle

# -------------------------------------------------------------------
# EOF analysis

def get_EOF_model(
    data_input: xr.DataArray, 
    time_dim: str = 'time', 
    n_modes = 10, 
    standardize = False, 
    use_coslat = False, 
    use_weights = False,
    weights = None,
    solver = 'auto',
    **solver_kwargs,
):
    """
    
    """
    from xeofs.models import EOF

    model = EOF(
        n_modes = n_modes, 
        standardize = standardize, 
        use_coslat = use_coslat, 
        use_weights = use_weights,
        solver = solver,
        solver_kwargs = solver_kwargs,
    )
    model.fit(data_input, dim = time_dim, weights = weights)

    return model

def calc_EOF_analysis(
    model: xeofs.models.eof.EOF,
    mini_summary = False,
):
    model_output = xr.Dataset()
    model_output['EOF'] = model.components()
    model_output['PC'] = model.scores()
    if mini_summary == True:
        return model_output
    model_output['explained_variance'] = model.explained_variance()
    model_output['explained_variance_ratio'] = model.explained_variance_ratio()
    model_output['singular_values'] = model.singular_values()

    return model_output

def get_EOF_projection(
    model: xeofs.models.eof.EOF,
    data: xr.DataArray,
):
    return model.transform(data)

def save_EOF_model(
    model: xeofs.models.eof.EOF,
    path: str,
):
    """
    
    """
    output_hal = open(path, 'wb')
    str = pickle.dumps(model)
    output_hal.write(str)
    output_hal.close()

def load_EOF_model(
    path: str,
):
    """
    
    """
    with open(path,'rb') as file:
        model = pickle.loads(file.read())
    return model

# -------------------------------------------------------------------
# MCA analysis

def get_MCA_model(
    data_input1: xr.DataArray, 
    data_input2: xr.DataArray, 
    time_dim: str = 'time', 
    n_modes = 10, 
    standardize = False, 
    use_coslat = False, 
    use_weights = False,
    weights1 = None,
    weights2 = None,
    n_pca_modes = None,
    solver = 'auto',
    **solver_kwargs
):
    """
    
    """
    from xeofs.models import MCA

    if n_pca_modes == None:
        n_pca_modes = data_input1[time_dim].shape[0]

    model = MCA(
        n_modes = n_modes, 
        standardize = standardize, 
        use_coslat = use_coslat, 
        use_weights = use_weights,
        solver = solver,
        n_pca_modes = n_pca_modes,
        solver_kwargs = solver_kwargs,
    )
    model.fit(data_input1, data_input2, dim = time_dim, weights1 = weights1, weights2 = weights2)

    return model

def calc_MCA_analysis(
    model: xeofs.models.mca.MCA,
    correction = None, 
    alpha = 0.05,
    mini_summary = False
):
    from datatree import DataTree

    model_output = DataTree(name="root")

    # components
    left_components = model.components()[0]
    right_components = model.components()[1]

    model_output['EOF/left_EOF'] = DataTree(left_components)
    model_output['EOF/right_EOF'] = DataTree(right_components)

    # scores
    left_scores = model.components()[0]
    right_scores = model.components()[1]

    model_output['PC/left_PC'] = DataTree(left_scores)
    model_output['PC/right_PC'] = DataTree(right_scores)

    if mini_summary == True:
        return model_output

    # covariance_fraction
    covariance_fraction = model.covariance_fraction()
    model_output['covariance_fraction'] = DataTree(covariance_fraction)

    # singular_values
    singular_values = model.singular_values()
    model_output['singular_values'] = DataTree(singular_values)

    # squared_covariance
    squared_covariance = model.squared_covariance()
    model_output['squared_covariance'] = DataTree(squared_covariance)

    # squared_covariance_fraction
    squared_covariance_fraction = model.squared_covariance_fraction()
    model_output['squared_covariance_fraction'] = DataTree(squared_covariance_fraction)

    # heterogeneous_patterns
    tmp = model.heterogeneous_patterns(correction = correction, alpha = alpha)
    left_heterogeneous_patterns = tmp[0][0]
    right_heterogeneous_patterns = tmp[0][1]
    pvalues_of_left_heterogeneous_patterns = tmp[1][0]
    pvalues_of_right_heterogeneous_patterns = tmp[1][1]

    model_output['heterogeneous_patterns/left_heterogeneous_patterns'] = left_heterogeneous_patterns
    model_output['heterogeneous_patterns/right_heterogeneous_patterns'] = right_heterogeneous_patterns
    model_output['heterogeneous_patterns/pvalues_of_left_heterogeneous_patterns'] = pvalues_of_left_heterogeneous_patterns
    model_output['heterogeneous_patterns/pvalues_of_right_heterogeneous_patterns'] = pvalues_of_right_heterogeneous_patterns

    # homogeneous_patterns
    tmp = model.homogeneous_patterns(correction = correction, alpha = alpha)
    left_homogeneous_patterns = tmp[0][0]
    right_homogeneous_patterns = tmp[0][1]
    pvalues_of_left_homogeneous_patterns = tmp[1][0]
    pvalues_of_right_homogeneous_patterns = tmp[1][1]

    model_output['homogeneous_patterns/left_homogeneous_patterns'] = left_homogeneous_patterns
    model_output['homogeneous_patterns/right_homogeneous_patterns'] = right_homogeneous_patterns
    model_output['homogeneous_patterns/pvalues_of_left_homogeneous_patterns'] = pvalues_of_left_homogeneous_patterns
    model_output['homogeneous_patterns/pvalues_of_right_homogeneous_patterns'] = pvalues_of_right_homogeneous_patterns

    return model_output

def get_MCA_projection(
    model: xeofs.models.mca.MCA,
    **kwargs,
):
    """
    
    """
    return model.transform(kwargs)

def save_MCA_model(
    model: xeofs.models.mca.MCA,
    path: str,
):
    """
    
    """
    output_hal = open(path, 'wb')
    str = pickle.dumps(model)
    output_hal.write(str)
    output_hal.close()

def load_MCA_model(
    path: str,
):
    """
    
    """
    with open(path,'rb') as file:
        model = pickle.loads(file.read())
    return model