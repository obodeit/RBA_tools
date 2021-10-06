# python 2/3 compatibility
from __future__ import division, print_function

# package imports
from rbatools.constraint_Info_Block import ConstraintBlock


class DensityConstraints(ConstraintBlock):
    """
    Class holding information on the constraints regarding the compartment-densities in the model.

   Attributes
   ----------
   Elements : dict
       Each model metabolite-constraint is represented by a key.
       The values, holding information on each process-constraint, are dicts with predefined keys:
           'ID' : metabolite ID in model (type str)
           'AssociatedCompartment' : ID of compartment this constraint defines the capacity for (type str)
           'Type' : Equality or inequality (type dict)
            'CapacityParameter' : Terms/Parameters defining efficiency (type list)
                Keys:
                    Aggregate parameters (efficiency is their product)
                Values:
                    'FunctionID': ID of parameter
                    'FunctionType': Type of mathematical function
                    'IndependentVariable' : variable their values depend on
                    'FunctionParameters' : parameters for function (depend on function type)
    """

    def fromFiles(self, model, Cs, matrix):
        index = 0
        self.Elements = {}
        for i in Cs['DensityConsts'].keys():
            index += 1
            if matrix.row_signs[Cs['DensityConsts'][i]] == 'L':
                cSign = '<='
            if matrix.row_signs[Cs['DensityConsts'][i]] == 'E':
                cSign = '='
            compartment = i.split('_density')[0]
            sizePar = getSizeParameter(model, compartment)
            self.Elements[i] = {'ID': i,
                                'AssociatedCompartment': compartment,
                                'index': index,
                                'CapacityParameter': getParameterFunction(model, sizePar),
                                'Type': cSign}


def getSizeParameter(model, compartment):
    s = ''
    for c in model.density.__dict__['target_densities'].__dict__['_elements']:
        if str(c.compartment) is compartment:
            s = c.__dict__['upper_bound']
    return(s)


def getParameterFunction(model, param):
    par = []
    if param in list(model.parameters.functions.__dict__['_elements_by_id']):
        par = [{param: getElementaryFunctionInfo(model, param)}]
    if param in list(model.parameters.aggregates.__dict__['_elements_by_id']):
        x = model.parameters.aggregates.__dict__['_elements_by_id'][param].__dict__
        type = x['type']
        elementaryFunctions = [i.__dict__['function'] for i in x['function_references']._elements]
        par = []
        for ef in elementaryFunctions:
            par.append({ef: getElementaryFunctionInfo(model, ef)})
    return(par)


def getElementaryFunctionInfo(model, funct):
    x = model.parameters.functions.__dict__['_elements_by_id'][funct].__dict__
    type = x['type']
    var = x['variable']
    if type == 'constant':
        F = {'C': str(x['parameters'].__dict__['_elements_by_id']['CONSTANT'].__dict__['value'])}
    elif type == 'exponential':
        F = {'Lambda': str(x['parameters'].__dict__['_elements_by_id']['RATE'].__dict__['value'])}
    elif type == 'indicator':
        F = {}
        if 'X_MIN' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'xMin': str(x['parameters'].__dict__[
                     '_elements_by_id']['X_MIN'].__dict__['value'])})
        if 'X_MAX' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'xMax': str(x['parameters'].__dict__[
                     '_elements_by_id']['X_MAX'].__dict__['value'])})
    elif type == 'linear':
        F = {}
        if 'LINEAR_COEF' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'A': str(x['parameters'].__dict__['_elements_by_id']
                               ['LINEAR_COEF'].__dict__['value'])})
        if 'LINEAR_CONSTANT' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'C': str(x['parameters'].__dict__['_elements_by_id']
                               ['LINEAR_CONSTANT'].__dict__['value'])})
        if 'X_MIN' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'xMin': str(x['parameters'].__dict__[
                     '_elements_by_id']['X_MIN'].__dict__['value'])})
        if 'X_MAX' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'xMax': str(x['parameters'].__dict__[
                     '_elements_by_id']['X_MAX'].__dict__['value'])})
        if 'Y_MIN' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'yMin': str(x['parameters'].__dict__[
                     '_elements_by_id']['Y_MIN'].__dict__['value'])})
        if 'Y_MAX' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'yMax': str(x['parameters'].__dict__[
                     '_elements_by_id']['Y_MAX'].__dict__['value'])})
    elif type == 'michaelisMenten':
        F = {}
        if 'kmax' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'Vmax': str(x['parameters'].__dict__[
                     '_elements_by_id']['kmax'].__dict__['value'])})
        if 'Km' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'Km': str(x['parameters'].__dict__[
                     '_elements_by_id']['Km'].__dict__['value'])})
        if 'Y_MIN' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'yMin': str(x['parameters'].__dict__[
                     '_elements_by_id']['Y_MIN'].__dict__['value'])})
    elif type == 'competitiveInhibition':
        F = {}
        if 'kmax' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'Vmax': str(x['parameters'].__dict__[
                     '_elements_by_id']['kmax'].__dict__['value'])})
        if 'Km' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'Km': str(x['parameters'].__dict__[
                     '_elements_by_id']['Km'].__dict__['value'])})
        if 'Y_MIN' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'yMin': str(x['parameters'].__dict__[
                     '_elements_by_id']['Y_MIN'].__dict__['value'])})
        if 'Ki' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'Ki': str(x['parameters'].__dict__[
                     '_elements_by_id']['Ki'].__dict__['value'])})
        if 'I' in list(x['parameters'].__dict__['_elements_by_id'].keys()):
            F.update({'I': str(x['parameters'].__dict__['_elements_by_id']['I'].__dict__['value'])})
    elif type == 'inverse':
        F = {'C': str(x['parameters'].__dict__['_elements_by_id']['CONSTANT'].__dict__['value'])}
    else:
        F = {}
    return({'FunctionID': funct, 'FunctionType': type, 'IndependentVariable': var, 'FunctionParameters': F})
