# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import json
import copy
import pandas
import numpy
from rbatools.information_block import InformationBlock


class DataBlock(InformationBlock):
    """
    Class holding data from model-simulations.

    """

    def addEntries(self, Dict):
        for i in Dict.keys():
            self.Elements[i] = Dict[i]

    def fromDict(self, Dict):
        self.Elements = Dict

    def JSONize(self):
        Block = self.Elements
        block2 = copy.deepcopy(Block)
        for i in list(Block.keys()):
            if type(Block[i]) is dict:
                for j in list(Block[i].keys()):
                    block2[i][j] = json.dumps(Block[i][j], default=JSON_Int64_compensation)
            else:
                block2[i] = json.dumps(Block[i], default=JSON_Int64_compensation)
        return(block2)

    def toDataFrame(self, Col_list=None):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            if Col_list is None:
                fields = list(Block[list(Block.keys())[0]].keys())
            else:
                fields = Col_list
            TableOut = pandas.DataFrame(index=list(Block.keys()), columns=fields)
            for i in list(Block.keys()):
                for j in fields:
                    TableOut.loc[i, j] = Block[i][j]
            return TableOut
        else:
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

            TableOut = pandas.DataFrame(columns=fields)
            for i in list(Block.keys()):
                for j in fields:
                    entry = Block[i][j]
                    intString = None
                    if isinstance(entry, str):
                        if len(entry) > 0:
                            intString = entry.replace("'", "")
                    elif isinstance(entry, list):
                        if len(entry) > 0:
                            intString = json.dumps(entry, default=JSON_Int64_compensation)
                    elif isinstance(entry, dict):
                        if len(list(entry.keys())) > 0:
                            intString = json.dumps(entry, default=JSON_Int64_compensation)
                    elif isinstance(entry, set):
                        if len(list(entry)) > 0:
                            intString = json.dumps(entry, default=JSON_Int64_compensation)
                    else:
                        intString = json.dumps(entry, default=JSON_Int64_compensation)
                    TableOut.loc[i, j] = intString

            if len(list(NameList)) == len(list(TableOut)):
                TableOut.columns = list(NameList)
            return TableOut
        else:
            return pandas.DataFrame(columns=NameList)

    def toSBtab(self, table_id, table_type, document_name=None, table_name=None, document=None, unit=None, Col_list=None, NameList=None):
        from sbtab import SBtab
        DF = self.toDataFrame(Col_list=Col_list)
        DF.reset_index(drop=False, inplace=True)
        DF.rename(columns={'index': 'VariableID'}, inplace=True)
        return(SBtab.SBtabTable.from_data_frame(df=DF, table_id=table_id, table_type=table_type, document_name=document_name, table_name=table_name, document=document, unit=unit, sbtab_version='1.0'))


def JSON_Int64_compensation(o):
    if isinstance(o, numpy.int64):
        return int(o)
    raise TypeError
