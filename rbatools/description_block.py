# python 2/3 compatibility
from __future__ import division, print_function

import copy
import json
import pandas

# package imports
from rbatools.information_block import InformationBlock


class DescriptionBlock(InformationBlock):
    """
    Class holding general information on model.

    Author, metabolic reconstruction etc...

    Attributes
    ----------
    Elements : Dictionary information on model.


    """

    def __init__(self):
        self.Elements = {}

    def fromFiles(self, File):
        for i in File.index.tolist():
            self.Elements[i] = File.loc[i][1]

    def addEntries(self, Dict):
        for i in Dict.keys():
            self.Elements[i] = Dict[i]

    def JSONize(self):
        Block = self.Elements
        block2 = copy.deepcopy(Block)
        for i in Block.keys():
            block2[i] = json.dumps(Block[i])
        return(block2)

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
                if isinstance(Block[i], str):
                    value_entry = Block[i].replace("'", "")
                Val = json.dumps(value_entry, default=JSON_Int64_compensation)
                TableOut.loc[i, 'Measure'] = '"'+Var.replace('"', '')+'"'
                TableOut.loc[i, 'Value'] = '"'+Val.replace('"', '')+'"'
            if len(list(NameList)) == len(list(TableOut)):
                TableOut.columns = list(NameList)
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

    def toDataFrame_RunInfo(self):
        Block = self.Elements
        if len(list(Block.keys())) > 0:
            runs = list(Block[list(Block.keys())[0]].keys())
            fields = ['Property']+runs
            TableOut = pandas.DataFrame(columns=fields, index=list(Block.keys()))
            for i in list(Block.keys()):
                Var = json.dumps(i, default=JSON_Int64_compensation)
                if "'" in Var:
                    Var = Var.replace("'", "")
                TableOut.loc[i, 'Property'] = Var
                for j in runs:
                    Val = json.dumps(Block[i][j], default=JSON_Int64_compensation)
                    if "'" in Val:
                        Val = Val.replace("'", "")
                    TableOut.loc[i, j] = Val.replace('"', '')
            return TableOut
        if len(list(Block.keys())) == 0:
            return pandas.DataFrame()

    def toSBtab_RunInfo_forDoc(self, table_id, table_type, document_name=None, table_name=None, document=None, unit=None, *NameList):
        SBtab_Colnames = []
        if len(list(NameList)) > 0:
            if len(list(NameList[0])) > 0:
                SBtab_Colnames = list(NameList[0])
        from sbtab import SBtab
#          print(SBtab_Colnames)
#          DF=self.toDataFrame_RunInfo(SBtab_Colnames)
        DF = self.toDataFrame_RunInfo()
        return(SBtab.SBtabTable.from_data_frame(df=DF, table_id=table_id, table_type=table_type, document_name=document_name, table_name=table_name, document=document, unit=unit, sbtab_version='1.0'))

    def toSBtab_RunInfo(self, table_id, table_type, document_name=None, table_name=None, document=None, unit=None):
        from sbtab import SBtab
        DF = self.toDataFrame_RunInfo()
        SB = SBtab.SBtabTable.from_data_frame(df=DF, table_id=table_id, table_type=table_type,
                                              document_name=document_name, table_name=table_name, document=document, unit=unit, sbtab_version='1.0')
        SB.write(table_type+'.tsv')

#     def toSBtab(self,document_name, table_type, table_name,document, unit):
#          from sbtab import SBtab
#          DF=self.toDataFrame()
#          SB=SBtab.SBtabTable.from_data_frame(DF, document_name, table_type, table_name, document, unit, sbtab_version='1.0')
#          SB.write(table_type)
#          f=open(document_name+'.tsv','w')
#          f.write(SB.table_string)
#          f.close()


def JSON_Int64_compensation(o):
    import numpy
    if isinstance(o, numpy.int64):
        return int(o)
    raise TypeError
