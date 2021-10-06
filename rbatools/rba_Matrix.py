from __future__ import division, print_function
import numpy
import scipy


class RBA_Matrix(object):
    """
    Class holding RBA-Linear Problem.
    Required Object type for linear problems (matrices), to be used in this RBA-API.
    Should include all fields, required for a linear problem in the COBRA format, as attributes.

    Attributes
    ----------
    A : scipy.sparse.coo_matrix
        Lefthandside of constraint Matrix (aka Constraint Matrix)
    b : numpy.array
        Righthandside of Constraint Matrix
    row_signs : list
        Type of constraints ('E' for equality 'L' for lower-or-equal inequality) --> Ax=b or Ax<=b
    f : numyp.array
        Objective function linear get_coefficients
    LB : numpy.array
        Lower bounds of decision-variables
    UB : numpy.array
        Upper bounds of decision-variables
    row_names : list
        Names of constraints
    col_names : list
        Names of decision-variables
    rowIndicesMap : dict
        Dictionary mapping constraint names to their numeric index (generated automatically)
    colIndicesMap : dict
        Dictionary mapping variable names to their numeric index (generated automatically)

    Methods
    ----------
    __init__()
        Initiates with empty Problem. Attributes can be changed manually.

    loadMatrix(matrix)
        Imports information for attributes from input-objects

    mapIndices():
        Generates rowIndicesMap and colIndicesMap from attributes.

    AtoLiL()
        Transforms constraint matrix (LHS) to scipy.sparse.lil_matrix format.

    AtoCOO()
        Transforms constraint matrix (LHS) to scipy.sparse.coo_matrix format.

    tofloat64()
        Transforms all numerical attributes to number type numpy.float64.
    """

    def __init__(self):
        """
        Initiates with empty Problem. Attributes can be changed manually.
        """

        self.A = scipy.sparse.coo.coo_matrix(numpy.array([]))
        self.b = numpy.array([])
        self.f = numpy.array([])
        self.LB = numpy.array([])
        self.UB = numpy.array([])
        self.row_signs = []
        self.row_names = []
        self.col_names = []
        self.mapIndices()

        def AtoLiL(self):
            self.A = self.A.tolil()

        def AtoCOO(self):
            self.A = self.A.tocoo()

    def loadMatrix(self, matrix):
        """
        Imports information from compatible object.

        Imports information for attributes from input-objects (named matrix)
        with the respectively named fields, required.
        Required fields are:
        'A','b','row_signs','f','LB','UB','row_names' and 'col_names'.

        Parameters
        ----------
        matrix : any LP-object holding the required fields
            The matrix with elements to be added
        """

        if checkForAttributes(matrix, ['A', 'b', 'f', 'LB', 'UB', 'row_signs', 'row_names', 'col_names']):
            if type(matrix.A) is scipy.sparse.coo.coo_matrix:
                self.A = matrix.A.astype('float64')
            elif type(matrix.A) is numpy.ndarray:
                self.A = scipy.sparse.coo.coo_matrix(matrix.A.astype('float64'))
            elif type(matrix.A) is scipy.sparse.lil_matrix:
                self.A = scipy.sparse.coo.coo_matrix(matrix.A.astype('float64'))
            else:
                print('A must be of type coo_matrix or array')
                return
            if type(matrix.b) is list or type(matrix.b) is numpy.ndarray:
                self.b = matrix.b.astype('float64')
            else:
                print('b must be of list')
                return
            if type(matrix.f) is numpy.ndarray:
                self.f = matrix.f.astype('float64')
            else:
                print('f must be of type array')
                return
            if type(matrix.LB) is numpy.ndarray:
                self.LB = matrix.LB.astype('float64')
            else:
                print('LB must be of type array')
                return
            if type(matrix.UB) is numpy.ndarray:
                self.UB = matrix.UB.astype('float64')
            else:
                print('UB must be of type array')
                return
            if type(matrix.row_signs) is list:
                self.row_signs = matrix.row_signs
            else:
                print('row_signs must be of list')
                return
            if type(matrix.row_names) is list:
                self.row_names = matrix.row_names
            else:
                print('row_names must be of list')
                return
            if type(matrix.col_names) is list:
                self.col_names = matrix.col_names
            else:
                print('col_names must be of list')
                return
        else:
            print('Input does not have all necessary elements')
            return
        self.mapIndices()

    def scaleLHS(self, factor):
        self.A = self.A*factor

    def mapIndices(self):
        """
        Generates rowIndicesMap and colIndicesMap from attributes.

        rowIndicesMap = Dictionary: {'constraint_name':index,...}
        colIndicesMap = Dictionary: {'variable_name':index,...}
        """
        self.rowIndicesMap = dict(zip(self.row_names, list(range(len(self.row_names)))))
        self.colIndicesMap = dict(zip(self.col_names, list(range(len(self.col_names)))))

    def AtoLiL(self):
        """
        Transforms constraint matrix (LHS) to scipy.sparse.lil_matrix format.
        """
        self.A = scipy.sparse.lil_matrix(self.A)

    def AtoCOO(self):
        """
        Transforms constraint matrix (LHS) to scipy.sparse.coo_matrix format.
        """
        self.A = scipy.sparse.coo_matrix(self.A)

    def tofloat64(self):
        """
        Transforms all numerical attributes to number type numpy.float64.
        """
        self.A = self.A.astype('float64')
        self.b = self.b.astype('float64')
        self.f = self.f.astype('float64')
        self.LB = self.LB.astype('float64')
        self.UB = self.UB.astype('float64')


def checkForAttributes(obj, attr):
    for i in attr:
        x = getattr(obj, i, None)
        if x is None:
            return(True)
    return(True)
