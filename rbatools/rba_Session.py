# python 2/3 compatibility
from __future__ import division, print_function
import sys
import os.path
import numpy
import pandas
import copy
import difflib
import scipy
import collections
import json
# package imports
import rba
from rbatools.rba_SimulationData import RBA_SimulationData
from rbatools.rba_SimulationParameters import RBA_SimulationParameters
from rbatools.rba_ModelStructure import RBA_ModelStructure
from rbatools.rba_Problem import RBA_Problem
from rbatools.rba_Matrix import RBA_Matrix
from rbatools.rba_LP import RBA_LP
from rbatools.rba_LogBook import RBA_LogBook


class RBA_Session(object):
    """
    Top level of the RBA API.

    Attributes
    ----------
    xml_dir : str
        Current Growth rate as numeric value
    model : rba.RbaModel
        Current Growth rate as numeric value
    matrices : rba.ConstraintMatrix
        Current Growth rate as numeric value
    solver : rba.Solver
        Current Growth rate as numeric value
    Problem : rbatools.RBA_Problem
        Current Growth rate as numeric value
    Medium : dict
        Current Growth rate as numeric value
    ModelStructure : rbatools.RBA_ModelStructure
        Current Growth rate as numeric value
    Results : dict
        Current Growth rate as numeric value
    Parameters : dict
        Current Growth rate as numeric value
    SimulationData : rbatools.RBA_SimulationData
        Current Growth rate as numeric value
    SimulationParameters : rbatools.RBA_SimulationParameters
        Current Growth rate as numeric value

    Methods
    ----------
    __init__(xml_dir)
        Creates controler object from files

    reloadModel()
        Reloads model from files

    recordResults(runName)
        Records Simulation output for further use.

    recordParameters(runName)
        Records Simulation parameters for further use.

    clearResults()
        Removes all previosly recorded results.

    clearParameters()
        Removes all previosly recorded parameters.

    writeResults(session_name='', digits=10)
        Creates SimulationData and SimulationParameters objects from recordings.

    setMu(Mu)
        Sets growth-rate to desired value.

    doSolve(runName='DontSave')
        Solves problem to find solution.

    findMaxGrowthRate(precision=0.00001, max=4, recording=False)
        Applies dichotomy-search to find the maximal feasible growth-rate.

    findMinMediumConcentration(metabolite, precision=0.00001, max=100, recording=False)
        Applies dichotomy-search to find the minimal feasible concentration of
        growth-substrate in medium.

    setMedium(changes)
        Sets the concentration of growth-substrate in medium.

    knockOut(gene)
        Emulates a gene knock out.

    FeasibleRange(variables)
        Determines the feasible range of model variables.

    ConstraintSaturation(constraints)
        Determines the saturation of model constraints at current solution.

    ParetoFront(variables, N)
        Determine Pareto front of two model variables.

    addProtein(input)
        Adds representation of individual proteins to problem.

    returnExchangeFluxes()

    """

    def __init__(self, xml_dir):
        """
        Creates RBA_Session object from files

        Parameters
        ----------
        xml_dir : str
            Path to the directory where rba-model files are located.
        """
        self.xml_dir = xml_dir
        self.LogBook = RBA_LogBook('Controler')

        if not hasattr(self, 'ModelStructure'):
            if os.path.isfile(str(self.xml_dir+'/ModelStructure.json')):
                self.ModelStructure = RBA_ModelStructure()
                with open(str(self.xml_dir+'/ModelStructure.json'), 'r') as myfile:
                    data = myfile.read()
                self.ModelStructure.fromJSON(inputString=data)
            else:
                self.build_ModelStructure()

        self.model = rba.RbaModel.from_xml(input_dir=xml_dir)
        self.matrices = rba.ConstraintMatrix(model=self.model)
        self.solver = rba.Solver(matrix=self.matrices)

        self.LogBook.addEntry('Model loaded from {}.'.format(self.xml_dir))
        self.Problem = RBA_Problem(solver=self.solver)

        medium = pandas.read_csv(xml_dir+'/medium.tsv', sep='\t')
        self.Medium = dict(zip(list(medium.iloc[:, 0]), [float(i)
                                                         for i in list(medium.iloc[:, 1])]))

        self.Mu = self.Problem.Mu
        self.ExchangeMap = buildExchangeMap(self)

    def build_ModelStructure(self):
        """
        Rebuilds model structure object from model files and stores as json.
        """
        self.ModelStructure = RBA_ModelStructure()
        self.ModelStructure.fromFiles(xml_dir=self.xml_dir)
        self.ModelStructure.exportJSON(path=self.xml_dir)

    def rebuild_from_model(self):
        """
        Rebuilds computational model-representation from own attribute "model" (rba.RbaModel-object).
        """
        self.LogBook.addEntry('Model rebuilt.')
        self.matrices = rba.ConstraintMatrix(model=self.model)
        self.solver = rba.Solver(matrix=self.matrices)
        self.Problem = RBA_Problem(solver=self.solver)
        self.setMedium(changes=self.Medium)

    def reloadModel(self):
        """
        Reloads model from xml-files and then rebuild computational model-representation.
        """
        self.LogBook.addEntry('Model reloaded from {}.'.format(self.xml_dir))
        self.model = rba.RbaModel.from_xml(input_dir=self.xml_dir)
        self.rebuild_from_model()

    def recordResults(self, runName):
        """
        Records Simulation output for further use.
        and strores them in own 'Results'-attribute as pandas.DataFrames in a dictionary with the respective run-name being a column in all DataFrames.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        self.LogBook.addEntry('Solution recorded under {}.'.format(runName))
        if not hasattr(self, 'Results'):
            self.Results = {'Reactions': pandas.DataFrame(index=list(self.ModelStructure.ReactionInfo.Elements.keys())),
                            'Enzymes': pandas.DataFrame(index=list(self.ModelStructure.EnzymeInfo.Elements.keys())),
                            'Processes': pandas.DataFrame(index=[self.ModelStructure.ProcessInfo.Elements[i]['ID']+'_machinery' for i in self.ModelStructure.ProcessInfo.Elements.keys()]),
                            'Proteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinMatrix['Proteins'])),
                            'ProtoProteins': pandas.DataFrame(index=list(self.ModelStructure.ProteinGeneMatrix['ProtoProteins'])),
                            'Constraints': pandas.DataFrame(index=self.Problem.LP.row_names),
                            'SolutionType': pandas.DataFrame(index=['SolutionType']),
                            'ObjectiveFunction': pandas.DataFrame(index=self.Problem.LP.col_names),
                            'Mu': pandas.DataFrame(index=['Mu']),
                            'ObjectiveValue': pandas.DataFrame(index=['ObjectiveValue']),
                            'ExchangeFluxes': pandas.DataFrame(index=list(self.ExchangeMap.keys()))}

        Exchanges = self.returnExchangeFluxes()
        for i in Exchanges.keys():
            self.Results['ExchangeFluxes'].loc[i, runName] = checkwithsolutionfeasibility(
                Value=Exchanges[i], Session=self)

        self.Results['Reactions'][runName] = [checkwithsolutionfeasibility(
            Value=self.Problem.SolutionValues[i], Session=self) for i in list(self.Results['Reactions'].index)]
        self.Results['Enzymes'][runName] = [checkwithsolutionfeasibility(
            Value=self.Problem.SolutionValues[i], Session=self) for i in list(self.Results['Enzymes'].index)]
        self.Results['Processes'][runName] = [checkwithsolutionfeasibility(
            Value=self.Problem.SolutionValues[i], Session=self) for i in list(self.Results['Processes'].index)]
        self.Results['Constraints'][runName] = [checkwithsolutionfeasibility(
            Value=self.Problem.DualValues[i], Session=self) for i in self.Problem.LP.row_names]
        self.Results['Proteins'][runName] = ProteomeRecording(self, runName)
        self.Results['ProtoProteins'][runName] = ProtoProteomeRecording(
            self, runName, self.Results['Proteins'])
        self.Results['SolutionType'][runName] = checkwithsolutionfeasibility(Value=self.Problem.SolutionType, Session=self)
        self.Results['Mu'][runName] = checkwithsolutionfeasibility(Value=self.Problem.Mu, Session=self)
        self.Results['ObjectiveValue'][runName] = checkwithsolutionfeasibility(
            Value=self.Problem.ObjectiveValue, Session=self)
        self.Results['ObjectiveFunction'][runName] = list(self.Problem.getObjective().values())

    def recordParameters(self, runName):
        """
        Records Simulation parameters (LP-coefficients etc.) for further use.
        and strores them in own 'Parameters'-attribute as pandas.DataFrames in a dictionary with the respective run-name being a column in all DataFrames.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all Data, originating from these.
        """
        self.LogBook.addEntry('Coefficients recorded under {}.'.format(runName))
        EnzymeCapacities = self.Problem.getEnzymeCapacities()
        ProcessCapacities = self.Problem.getProcessCapacities()
        CompartmentCapacities = self.Problem.getCompartmentCapacities()
        if not hasattr(self, 'Parameters'):
            self.Parameters = {'EnzymeEfficiencies_FW': pandas.DataFrame(index=list(EnzymeCapacities.keys())),
                               'EnzymeEfficiencies_BW': pandas.DataFrame(index=list(EnzymeCapacities.keys())),
                               'NetProcessEfficiencies': pandas.DataFrame(index=list(ProcessCapacities.keys())),
                               'CompartmentCapacities': pandas.DataFrame(index=list(CompartmentCapacities.keys()))}

        self.Parameters['EnzymeEfficiencies_FW'][runName] = [
            EnzymeCapacities[i]['Forward'] for i in list(EnzymeCapacities.keys())]
        self.Parameters['EnzymeEfficiencies_BW'][runName] = [
            EnzymeCapacities[i]['Backward'] for i in list(EnzymeCapacities.keys())]
        self.Parameters['NetProcessEfficiencies'][runName] = [ProcessCapacities[i]
                                                              for i in list(ProcessCapacities.keys())]
        self.Parameters['CompartmentCapacities'][runName] = [CompartmentCapacities[i]
                                                             for i in list(CompartmentCapacities.keys())]

    def clearResults(self):
        """
        Removes all previosly recorded results and deletes own 'Results'-attribute.
        """
        self.LogBook.addEntry('Results cleared.')

        delattr(self, 'Results')

    def clearParameters(self):
        """
        Removes all previosly recorded parameters and deletes own 'Parameters'-attribute.
        """

        self.LogBook.addEntry('Parameters cleared.')
        delattr(self, 'Parameters')

    def writeResults(self, session_name='', digits=10, loggingIntermediateSteps=False):
        """
        Creates SimulationData and SimulationParameters objects from recordings ('Results'.'Parameters').

        Stores them as rbatools.RBA_SimulationData
        and rbatools.RBA_SimulationParameters objects as attributes.
        Access via attributes .SimulationData and SimulationParameters respectively.

        Parameters
        ----------
        digits : int
            Number of decimal places in the numeric results
            Default: 10
        session_name : str
            Name of Simulation session.
            Default: ''
        """
        self.LogBook.addEntry('Data written under {}.'.format(session_name))
        if hasattr(self, 'Results'):
            self.Results['uniqueReactions'] = mapIsoReactions(Controller=self)
            self.Results['Mu'] = self.Results['Mu'].round(digits)
            self.Results['ObjectiveValue'] = self.Results['ObjectiveValue'].round(digits)
            self.Results['Proteins'] = self.Results['Proteins'].round(digits)
            self.Results['uniqueReactions'] = self.Results['uniqueReactions'].round(digits)
            self.Results['Reactions'] = self.Results['Reactions'].round(digits)
            self.Results['Enzymes'] = self.Results['Enzymes'].round(digits)
            self.Results['Processes'] = self.Results['Processes'].round(digits)
            self.Results['Constraints'] = self.Results['Constraints'].round(digits)
            self.Results['ExchangeFluxes'] = self.Results['ExchangeFluxes'].round(digits)

            self.SimulationData = RBA_SimulationData(StaticData=self.ModelStructure)
            self.SimulationData.fromSimulationResults(Controller=self, session_name=session_name)

        if hasattr(self, 'Parameters'):
            self.Parameters['EnzymeEfficiencies_FW'] = self.Parameters['EnzymeEfficiencies_FW'].round(
                digits)
            self.Parameters['EnzymeEfficiencies_BW'] = self.Parameters['EnzymeEfficiencies_BW'].round(
                digits)
            self.Parameters['NetProcessEfficiencies'] = self.Parameters['NetProcessEfficiencies'].round(
                digits)
            self.Parameters['CompartmentCapacities'] = self.Parameters['CompartmentCapacities'].round(
                digits)
            self.SimulationParameters = RBA_SimulationParameters(StaticData=self.ModelStructure)
            self.SimulationParameters.fromSimulationResults(Controller=self)

    def returnExchangeFluxes(self):
        """
        Generates a dictonary with the exchang-rates of boundary-metabolites.

        Returns
        -------
        Dictonary with exchange-keys and respective -rates.
        """
        out = {}
        for j in self.ExchangeMap.keys():
            netflux = 0
            for k in self.ExchangeMap[j].keys():
                netflux += self.ExchangeMap[j][k]*self.Problem.SolutionValues[k]
            if netflux != 0:
                out[j] = netflux
        return(out)

    def setMu(self, Mu, loggingIntermediateSteps=False):
        """
        Sets growth-rate to desired value.

        Parameters
        ----------
        Mu : float
            Growth rate
        """
        self.LogBook.addEntry('Growth-rate changed:{} --> {}'.format(self.Mu, float(Mu)))
        self.Problem.setMu(Mu=float(Mu), ModelStructure=self.ModelStructure,
                           logging=loggingIntermediateSteps)
        self.Mu = float(Mu)

    def doSolve(self, runName='DontSave', feasibleStatuses=[1], try_unscaling_if_sol_status_is_five=True, loggingIntermediateSteps=False):
        """
        Solves problem to find solution.

        Does the same as rbatools.RBA_Problem.solveLP().
        Just has some automatic option for results-recording.

        Parameters
        ----------
        runName : str
            Name of observation.
            Serves as ID for all data, originating from this run.
            Special values :
                'DontSave' : Results are not recorded
                'Auto' : Results are automatically recorded
                         and appended to existing ones.
                    Named with number.
                Any other string: Results are recorded under this name.
            Default: 'DontSave'
        feasibleStatuses : list of int
            List with identifiers of acceptable solution statuses.
            (consult ILOG-CPLEX documentation for information on them).
            Default: [1]
        try_unscaling_if_sol_status_is_five : bool
        	If true; the problem will be attempted to be solved without scaling, 
        	if the scaled problem is feasible but the solution is not feasible 
        	after unscaling (CPLEX solution-status 5).
        	Default: True    
        """

        self.Problem.solveLP(feasibleStatuses=feasibleStatuses, try_unscaling_if_sol_status_is_five=try_unscaling_if_sol_status_is_five, logging=loggingIntermediateSteps)
        if self.Problem.Solved:
            if runName is not 'DontSave':
                if runName is 'Auto':
                    if hasattr(self, 'Results'):
                        name = str(self.Results['Reactions'].shape[1]+1)
                    if not hasattr(self, 'Results'):
                        name = '1'
                if runName is not 'Auto':
                    name = runName
                self.recordResults(runName=name)

    def findMaxGrowthRate(self, precision=0.001, max=4, start_value=None, recording=False, loggingIntermediateSteps=False, omit_objective=False, feasibleStatuses=[1], try_unscaling_if_sol_status_is_five=True):
        """
        Applies dichotomy-search to find the maximal feasible growth-rate.

        Parameters
        ----------
        precision : float
            Numberic precision with which maximum is approximated.
            Default : 0.00001
        max : float
            Defines the highest growth rate to be screened for.
            Default=4
        recording : bool
            Records intermediate feasible solutions
            while approaching the maximum growth-rate.
            Default : False
        feasibleStatuses : list of int
            List with identifiers of acceptable solution statuses.
            (consult ILOG-CPLEX documentation for information on them).
            Default: feasibleStatuses=[1]
        try_unscaling_if_sol_status_is_five : bool
        	If true; the problem will be attempted to be solved without scaling, 
        	if the scaled problem is feasible but the solution is not feasible 
        	after unscaling (CPLEX solution-status 5).
        	Default: try_unscaling_if_sol_status_is_five=True    

        Returns
        -------
        maximum feasible growth rate as float.
        """

        minMu = 0
        maxMu = max
        if start_value is None:
            testMu = maxMu
        else:
            testMu = start_value
        iteration = 0

        if omit_objective:
            old_Obj = self.Problem.getObjective()
            self.Problem.clearObjective()

        while (maxMu - minMu) > precision:
            self.setMu(Mu=testMu)
            self.Problem.solveLP(feasibleStatuses=feasibleStatuses,try_unscaling_if_sol_status_is_five=try_unscaling_if_sol_status_is_five,logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                iteration += 1
                if recording:
                    self.recordResults('DichotomyMu_iteration_'+str(iteration))
                minMu = testMu
            else:
                maxMu = testMu
            testMu = numpy.mean([maxMu, minMu])
        self.LogBook.addEntry('Maximal growth-rate found to be: {}.'.format(minMu))
        if minMu == max:
            print('Warning: Maximum growth rate might exceed specified range. Try rerunning this method with larger max-argument.')

        if omit_objective:
            self.Problem.setObjectiveCoefficients(old_Obj)
        self.setMu(Mu=minMu)
        self.Problem.solveLP(feasibleStatuses=feasibleStatuses,try_unscaling_if_sol_status_is_five=try_unscaling_if_sol_status_is_five,logging=False)
        return(minMu)

    def setMedium(self, changes, loggingIntermediateSteps=False):
        """
        Sets the concentration of specified growth-substrate(s) in medium.

        Parameters
        ----------
        changes : dict
            Keys : ID of metabolite(s) in medium.
            Values : New concention(s)
        """

        for species in (changes.keys()):
            self.Medium[species] = float(changes[species])

        self.Problem.ClassicRBAmatrix.set_medium(self.Medium)
        self.Problem.ClassicRBAmatrix.build_matrices(self.Mu)

        inputMatrix = RBA_Matrix()
        inputMatrix.loadMatrix(matrix=self.Problem.ClassicRBAmatrix)
        self.Problem.LP.updateMatrix(matrix=inputMatrix, Ainds=MediumDependentCoefficients_A(
            self), Binds=[], CTinds=[], LBinds=None, UBinds=None)

    def FeasibleRange(self, variables=None, loggingIntermediateSteps=False):
        """
        Determines the feasible range of model variables.

        Parameters
        ----------
        variables : str or list of str
            Specifies variable(s) for which the feasible range is to be determined.
            Optional input:
                If not provided all model-variables are taken

        Returns
        -------
        Dictionary with variable-names as keys and other dictionaries as values.
        The 'inner' dictionaries hold keys 'Min' and 'Max'
        with values representing lower and upper bound of feasible range respectively.
        E.g. : {'variableA':{'Min':42 , 'Max':9000},
                'variableB':{'Min':-9000 , 'Max':-42}}
        """

        if variables is not None:
        	if isinstance(variables, list):
        		VariablesInQuestion=variables
        	elif isinstance(variables, str):
        		VariablesInQuestion=[variables]   
        else:
        	VariablesInQuestion = self.Problem.LP.col_names
        			
        out = {}
        for i in VariablesInQuestion:
            min = numpy.nan
            max = numpy.nan
            self.Problem.clearObjective(logging=loggingIntermediateSteps)
#            self.Problem.setObjectiveCoefficients(inputDict=dict(
#                zip(self.Problem.LP.col_names, [0.0]*len(self.Problem.LP.col_names))))
            self.Problem.setObjectiveCoefficients(
                inputDict={i: 1.0}, logging=loggingIntermediateSteps)
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                min = self.Problem.SolutionValues[i]
            self.Problem.setObjectiveCoefficients(
                inputDict={i: -1.0}, logging=loggingIntermediateSteps)
            self.Problem.solveLP(logging=loggingIntermediateSteps)
            if self.Problem.Solved:
                max = self.Problem.SolutionValues[i]
            out.update({i: {'Min': min, 'Max': max}})
            self.LogBook.addEntry(
                'Feasible-range of {} determined to be between {} and {}.'.format(i, min, max))

        return(out)

    def ConstraintSaturation(self, constraints=None):
        """
        Determines the saturation of model constraints at current solution.

        Parameters
        ----------
        constraints : str or list of str
            Specifies constraints(s) for which the saturation is to be determined.
            Optional input:
                If not provided all model-constraints are taken

        Returns
        -------
        Pandas DataFrame with constraint-names as indices and the columns 'LHS', 'RHS', and 'Saturation'.
        	'LHS': The sum over the respoctive constraint-row multiplied elementwise with the solution vector.
        	'RHS': The value of the problem's righthand side, correesponding to the respective constraint.
        	'Saturation': The saturation of the respective constraint ('LHS'/'RHS').
        """
        if constraints is not None:
        	if isinstance(constraints, list):
        		ConstraintsInQuestion=constraints
        	elif isinstance(constraints, str):
        		ConstraintsInQuestion=[constraints]   
        else:
        	ConstraintsInQuestion = self.Problem.LP.row_names

        rhs = self.Problem.getRighthandSideValue(ConstraintsInQuestion)
        lhs = self.Problem.calculateLefthandSideValue(ConstraintsInQuestion)
        Out = pandas.DataFrame(columns=['LHS', 'RHS', 'Saturation'], index=ConstraintsInQuestion)
        for i in ConstraintsInQuestion:
            lhval = lhs[i]
            rhval = rhs[i]
            sat = numpy.nan
            if rhval != 0:
                sat = lhval/rhval
            Out.loc[i, 'LHS'] = lhval
            Out.loc[i, 'RHS'] = rhval
            Out.loc[i, 'Saturation'] = sat
            self.LogBook.addEntry(
                'Saturation of constraint {} determined to be {}.'.format(i, sat))
        return(Out)

    def addExchangeReactions(self):
        """
        Adds explicit exchange-reactions of boundary-metabolites to RBA-problem, 
        named R_EX_ followed by metabolite name (without M_ prefix).
        """
        Mets_external = [m.id for m in self.model.metabolism.species if m.boundary_condition]
        Mets_internal = [m.id for m in self.model.metabolism.species if not m.boundary_condition]
        Reactions = [r.id for r in self.model.metabolism.reactions]
        full_S = rba.core.metabolism.build_S(
            Mets_external+Mets_internal, self.model.metabolism.reactions)
        S_M_ext = full_S[:len(Mets_external), ].toarray()
        col_indices_toremove = []
        for i in range(S_M_ext.shape[1]):
            s_col_uniques = list(set(list(S_M_ext[:, i])))
            if len(s_col_uniques) == 1:
                if s_col_uniques[0] == 0:
                    col_indices_toremove.append(i)
        RemainingReactions = [i for i in Reactions if Reactions.index(
            i) not in col_indices_toremove]
        S_ext = numpy.delete(S_M_ext, col_indices_toremove, axis=1)
        A = numpy.concatenate((S_ext, numpy.eye(len(Mets_external))), axis=1, out=None)
        ColNames = RemainingReactions+[str('R_EX_'+i.split('M_')[-1]) for i in Mets_external]
        # print(str('R_EX_'+i.split('M_')[-1]))
        LBs = list([self.Problem.LP.LB[self.Problem.LP.col_names.index(i)]
                    for i in RemainingReactions]+[-10000]*len(Mets_external))
        UBs = list([self.Problem.LP.UB[self.Problem.LP.col_names.index(i)]
                    for i in RemainingReactions]+[10000]*len(Mets_external))
        b = [0]*len(Mets_external)
        f = list([self.Problem.LP.f[self.Problem.LP.col_names.index(i)]
                  for i in RemainingReactions]+[0]*len(Mets_external))

        ExchangeMatrix = RBA_Matrix()
        ExchangeMatrix.A = scipy.sparse.coo_matrix(A)
        ExchangeMatrix.b = numpy.array([0]*len(Mets_external))
        ExchangeMatrix.f = numpy.array(f)
        ExchangeMatrix.LB = numpy.array(LBs)
        ExchangeMatrix.UB = numpy.array(UBs)
        ExchangeMatrix.row_signs = ['E']*len(Mets_external)
        ExchangeMatrix.row_names = Mets_external
        ExchangeMatrix.col_names = ColNames
        ExchangeMatrix.mapIndices()
        self.Problem.LP.addMatrix(matrix=ExchangeMatrix)

        self.ExchangeReactionMap = dict(
            zip(Mets_external, [str('R_EX_'+i.split('M_')[-1]) for i in Mets_external]))


def MediumDependentCoefficients_A(Controler):
    out = {}
    MedDepRxns = [list(i.keys()) for i in list(Controler.ExchangeMap.values())]
    MedDepRxnsFlatted = list(set([item for sublist in MedDepRxns for item in sublist]))
    for i in Controler.ModelStructure.EnzymeConstraintsInfo.Elements.keys():
        if Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['AssociatedReaction'] in MedDepRxnsFlatted:
            nonConst = False
            for j in Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['CapacityParameter']:
                if list(j.values())[0]['FunctionType'] != 'constant':
                    nonConst = True
            if nonConst:
                if Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['AssociatedReaction'] in list(out.keys()):
                    out[Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]
                        ['AssociatedReaction']].append(i)
                else:
                    out.update(
                        {Controler.ModelStructure.EnzymeConstraintsInfo.Elements[i]['AssociatedReaction']: [i]})
    return([(out[i][0], Controler.ModelStructure.ReactionInfo.Elements[i]['Enzyme'])for i in out.keys()])


def QualitativeMediumChange(Controller, changes, species):
    QualitativeMediumChange = False
    if float(Controller.Medium[species]) == float(0):
        if float(changes[species]) != float(0):
            boundValue = 1000.0
            QualitativeMediumChange = True
        else:
            return([QualitativeMediumChange])
    if float(Controller.Medium[species]) != float(0):
        if float(changes[species]) == float(0):
            boundValue = 0.0
            QualitativeMediumChange = True
        else:
            return([QualitativeMediumChange])
    return([QualitativeMediumChange, float(boundValue)])


def findExchangeReactions(Controller, species):
    Reactions = list(Controller.ExchangeMap[species].keys())
    exchanges = {}
    for i in Reactions:
        if Controller.ExchangeMap[species][i] > 0:
            exchanges.update({i: 'Product'})
        elif Controller.ExchangeMap[species][i] < 0:
            exchanges.update({i: 'Reactant'})
    return(exchanges)


def determineCoefficient(x, changes, species):
    multiplicativeFactors = []
    for k in x:
        result = 1
        type = list(k.values())[0]['FunctionType']
        pars = list(k.values())[0]['FunctionParameters']
        if type == 'constant':
            result = numpy.float64(pars['C'])
        if type == 'exponential':
            L = 1
            if 'Lambda' in list(pars.keys()):
                L = numpy.float64(pars['Lambda'])
            result = numpy.exp(float(changes[species])*L)
        if type == 'indicator':
            maxi = numpy.inf
            mini = -numpy.inf
            if 'xMax' in list(pars.keys()):
                maxi = numpy.float64(pars['xMax'])
            if 'xMin' in list(pars.keys()):
                mini = numpy.float64(pars['xMin'])
            result = (float(changes[species]) > mini) and (float(changes[species]) < maxi)
        if type == 'linear':
            X_maxi = numpy.inf
            X_mini = -numpy.inf
            Y_maxi = numpy.inf
            Y_mini = -numpy.inf
            A = 1
            C = 0
            if 'A' in list(pars.keys()):
                A = numpy.float64(pars['A'])
            if 'C' in list(pars.keys()):
                C = numpy.float64(pars['C'])
            if 'xMin' in list(pars.keys()):
                X_mini = numpy.float64(pars['xMin'])
            if 'xMax' in list(pars.keys()):
                X_maxi = numpy.float64(pars['xMax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            if 'yMax' in list(pars.keys()):
                Y_maxi = numpy.float64(pars['yMax'])
            X = float(changes[species])
            if float(changes[species]) < X_mini:
                X = X_mini
            if float(changes[species]) > X_maxi:
                X = X_maxi
            Y = A*X + C
            result = Y
            if Y < Y_mini:
                result = Y_mini
            if Y > Y_maxi:
                result = Y_maxi
        if type == 'michaelisMenten':
            Y_mini = -numpy.inf
            KM = 0
            VM = 1
            if 'Km' in list(pars.keys()):
                KM = numpy.float64(pars['Km'])
            if 'Vmax' in list(pars.keys()):
                VM = numpy.float64(pars['Vmax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            Y = VM*float(changes[species])/(float(changes[species])+KM)
            result = Y
            if Y < Y_mini:
                result = Y_mini
        if type == 'competitiveInhibition':
            Y_mini = -numpy.inf
            KM = 0
            VM = 1
            KI = 0
            I = 0
            if 'Ki' in list(pars.keys()):
                KI = numpy.float64(pars['Ki'])
            if 'I' in list(pars.keys()):
                I = numpy.float64(pars['I'])
            if 'Km' in list(pars.keys()):
                KM = numpy.float64(pars['Km'])
            if 'Vmax' in list(pars.keys()):
                VM = numpy.float64(pars['Vmax'])
            if 'yMin' in list(pars.keys()):
                Y_mini = numpy.float64(pars['yMin'])
            Y = VM*float(changes[species])/(float(changes[species])+KM*(1+I/KI))
            result = Y
            if Y < Y_mini:
                result = Y_mini
        if type == 'inverse':
            C = 1
            if 'C' in list(pars.keys()):
                C = numpy.float64(pars['C'])
            result = 1
            if float(changes[species]) is not 0:
                result = C/float(changes[i])
        multiplicativeFactors.append(result)
        value = numpy.prod(numpy.array(multiplicativeFactors))
    return(float(value))


def ProtoProteomeRecording(Controller, run, Proteinlevels):
    out = []
    for i in list(Controller.ModelStructure.ProteinGeneMatrix['ProtoProteins']):
        row_ind = list(Controller.ModelStructure.ProteinGeneMatrix['ProtoProteins']).index(i)
        # print(row_ind)
        nonZero = list(numpy.nonzero(
            Controller.ModelStructure.ProteinGeneMatrix['Matrix'][row_ind, :])[0])
        level = 0
        for j in nonZero:
            id = Controller.ModelStructure.ProteinGeneMatrix['Proteins'][j]
            level += Proteinlevels.loc[id, run]
        out.append(level)
    return(out)


def ProteomeRecording(Controller, run):

    EnzDF = pandas.DataFrame(index=Controller.Problem.Enzymes)
    PrcDF = pandas.DataFrame(index=Controller.Problem.Processes)
    EnzDF[run] = [Controller.Problem.SolutionValues[i]for i in Controller.Problem.Enzymes]
    PrcDF[run] = [Controller.Problem.SolutionValues[i]for i in Controller.Problem.Processes]

    ProteinProteinMatrix = numpy.array(
        Controller.ModelStructure.ProteinMatrix['Matrix']).astype(numpy.float64)
    C = Controller.ModelStructure.ProteinMatrix['Consumers']
    Consumers = []
    for i in C:
        if i.startswith('P_'):
            # Consumers.append(str(i+'_machinery'))
            Consumers.append(str(i))
        if not i.startswith('P_'):
            Consumers.append(i)
    Proteins = Controller.ModelStructure.ProteinMatrix['Proteins']
    DF = pandas.concat([EnzDF, PrcDF], axis=0)
    ProteinLevels = pandas.DataFrame(index=Proteins)
    vector = numpy.nan_to_num(DF[run].reindex(Consumers))
    Level = ProteinProteinMatrix.dot(vector)
    ProteinLevels[run] = Level
    addedProts = [col for col in Controller.Problem.LP.col_names if col.startswith('TotalLevel_')]
    if len(addedProts) > 0:
        for p in addedProts:
            protID = p.split('TotalLevel_')[1]
            ProteinLevels[run].loc[protID] = Controller.Problem.SolutionValues[p]
    return(list(ProteinLevels[run]))


def mapIsoReactions(Controller):
    if hasattr(Controller, 'Results'):
        out = pandas.DataFrame()
        for run in list(Controller.Results['Reactions'].columns):
            rf = dict(zip(list(Controller.Results['Reactions'].index), list(
                Controller.Results['Reactions'][run])))
            rf = {k: v for k, v in rf.items() if v != 0.}
            rf_merged = collections.defaultdict(float)
            for reac_id, flux_val in rf.items():
                if "duplicate" in reac_id:
                    last_idx = reac_id.index('duplicate') - 1
                    rf_merged[reac_id[:last_idx]] += flux_val
                else:
                    rf_merged[reac_id] += flux_val
            if len(list(out)) == 0:
                out[run] = list(rf_merged.values())
                out.index = list(rf_merged.keys())
            else:
                runDF = pandas.DataFrame(list(rf_merged.values()),
                                         index=list(rf_merged.keys()), columns=[run])
                runDF = runDF.reindex(list(set(list(out.index)).union(
                    set(list(rf_merged.keys())))), fill_value=0)
                out = out.reindex(list(set(list(out.index)).union(
                    set(list(rf_merged.keys())))), fill_value=0)
                out = out.join(runDF, how='outer')
        return(out)


def buildExchangeMap(Controller):
    """
    Returns a map of all metabolites, the corresponding transport-reactions and stoichiometires;
    exchanged with the medium.
    {Metabolite1 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2},
    {Metabolite2 : {ExchangeReaction1 : stoch-coefficient1 , ExchangeReaction2 : stoch-coefficient2}}

    Metabolite1 - ... MetaboliteN : All metabolite-species in the medium (see medium.tsv file)
    ExchangeReaction1 - ... ExchangeReactionN : All metabolic reactions, which exchange the respective metabolite with the medium.
    stoch-coefficient : Stochiometric coefficient with which the respective metabolite is exchanged by the corresponding reaction.
    (Negative when reaction transports metabolite out of the cell; and positive when inside the cell.)

    Parameters
    ----------
    Controller : rbatools.NewControler.RBA_newControler

    Returns
    -------
    Dict.
    """
    BoundaryMetabolites = [i for i in list(Controller.ModelStructure.MetaboliteInfo.Elements.keys(
    )) if Controller.ModelStructure.MetaboliteInfo.Elements[i]['boundary']]
    ExchangeMap = {}
    for bM in BoundaryMetabolites:
        for rxn in Controller.ModelStructure.MetaboliteInfo.Elements[bM]['ReactionsInvolvedWith']:
            Reactants = list(
                Controller.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'].keys())
            Products = list(Controller.ModelStructure.ReactionInfo.Elements[rxn]['Products'].keys())
            if len(list(set(list(Reactants+Products)))) > 1:
                for met in list(set(list(Reactants+Products))):
                    # if met != bM:
                    if met == bM:
                        MediumSpecies = findExchangeMetInMedium(met, Controller.Medium)
                        if met in Reactants:
                            stochCoeff = - \
                                Controller.ModelStructure.ReactionInfo.Elements[rxn]['Reactants'][met]
                        elif met in Products:
                            stochCoeff = Controller.ModelStructure.ReactionInfo.Elements[rxn]['Products'][met]
                        if MediumSpecies in list(ExchangeMap.keys()):
                            ExchangeMap[MediumSpecies].update({rxn: stochCoeff})
                        else:
                            ExchangeMap[MediumSpecies] = {rxn: stochCoeff}
    return(ExchangeMap)


def findExchangeMetInMedium(metabolite, Medium):
    """
    Returns the most likely species in the Medium, for any Metabolic species.
    Parameters
    ----------
    metabolite : str
    Medium : dict
    -------
    Most likely ID as str
    """
    if metabolite.endswith('_e'):
        out = difflib.get_close_matches('_e'.join(metabolite.split('_e')[:-1]), Medium, 1)
    else:
        out = difflib.get_close_matches(metabolite, Medium, 1)
    if len(out) > 0:
        return(out[0])
    else:
        return('')


def checkwithsolutionfeasibility(Value, Session):
    if Session.Problem.Solved:
        return(Value)
    else:
        return(numpy.nan)
