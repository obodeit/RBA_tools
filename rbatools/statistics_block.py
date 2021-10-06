# python 2/3 compatibility
from __future__ import division, print_function

import copy
import json
import pandas

from rbatools.information_block import InformationBlock


class StatisticsBlock(InformationBlock):
    """
    Class holding model-statistics.

    Brief summary of key-numbers of model.

    Attributes
    ----------
    Elements : Dictionary different numbers on model.

    """

    def derive(self, Dict):
        self.Elements = Dict

    def toDataFrame(self):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            fields = ['Measure', 'Value']
            TableOut = pandas.DataFrame(columns=fields, index=list(Block.keys()))
            for i in list(Block.keys()):
                Var = json.dumps(i, default=JSON_Int64_compensation)
                Val = json.dumps(Block[i], default=JSON_Int64_compensation)
                TableOut.loc[i, 'Measure'] = Var
                TableOut.loc[i, 'Value'] = Val
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

    def toDataFrame_SBtabCompatibility(self, NameList=None, Col_list=None):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            if Col_list is None:
                fields = list(Block[list(Block.keys())[0]].keys())
            else:
                fields = Col_list
            if NameList is not None:
                if len(fields) == len(NameList):
                    colNames = NameList
                else:
                    colNames = fields
            else:
                colNames = fields

            TableOut = pandas.DataFrame(columns=fields, index=list(Block.keys()))
            for i in list(Block.keys()):
                Var = json.dumps(i, default=JSON_Int64_compensation)
                Val = json.dumps(Block[i], default=JSON_Int64_compensation)
                if "'" in Var:
                    Var = Var.replace("'", "")
                if "'" in Val:
                    Val = Val.replace("'", "")
                TableOut.loc[i, 'Measure'] = Var
                TableOut.loc[i, 'Value'] = Val
# WORKS                    TableOut.loc[i,'Measure']='"'+Var+'"'
# WORKS                    TableOut.loc[i,'Value']='"'+Val+'"'
            if len(list(NameList)) == len(list(TableOut)):
                TableOut.columns = list(NameList)
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()


def JSON_Int64_compensation(o):
    if isinstance(o, numpy.int64):
        return int(o)
    raise TypeError
