# python 2/3 compatibility
from __future__ import division, print_function
import sys
import os.path
import numpy
import pandas
import copy
import json
import jxmlease
import xml.etree.ElementTree as ET
from sbtab import SBtab

# package imports
import rba
from rbatools.parameter_block import ParameterBlock


class RBA_SimulationParameters(object):
    """
    Class holding information on simulations with the model.

    Attributes
    ----------
    EnzymeCapacities : rbatools.ParameterBlock object holding enzyme capacities.

    ProcessCapacities : rbatools.ParameterBlock object holding process capacities.

    CompartmentCapacities : rbatools.ParameterBlock object holding compartment capacities.


    Methods
    ----------
    __init__:

    fromSimulationResults:
    Arguments: 'Controller'= rbatools.Controler object

    exportSBtab_OneFile

    """

    def __init__(self, StaticData):

        self.EnzymeCapacities_FW = ParameterBlock()
        self.EnzymeCapacities_BW = ParameterBlock()
        self.ProcessCapacities = ParameterBlock()
        self.CompartmentCapacities = ParameterBlock()
        self.EnzymeCapacities_FW.fromDict({})
        self.EnzymeCapacities_BW.fromDict({})
        self.ProcessCapacities.fromDict({})
        self.CompartmentCapacities.fromDict({})

    def fromSimulationResults(self, Controller):

        for effi in list(Controller.Parameters['EnzymeEfficiencies_FW'].index):
            if effi not in self.EnzymeCapacities_FW.Elements:
                self.EnzymeCapacities_FW.Elements.update({effi: {}})
            self.EnzymeCapacities_FW.Elements[effi].update({'ID': effi})
            for run in list(Controller.Parameters['EnzymeEfficiencies_FW']):
                self.EnzymeCapacities_FW.Elements[effi].update(
                    {run: json.dumps(Controller.Parameters['EnzymeEfficiencies_FW'].loc[effi, run])})

        for effi in list(Controller.Parameters['EnzymeEfficiencies_BW'].index):
            if effi not in self.EnzymeCapacities_BW.Elements:
                self.EnzymeCapacities_BW.Elements.update({effi: {}})
            self.EnzymeCapacities_BW.Elements[effi].update({'ID': effi})
            for run in list(Controller.Parameters['EnzymeEfficiencies_BW']):
                self.EnzymeCapacities_BW.Elements[effi].update(
                    {run: json.dumps(Controller.Parameters['EnzymeEfficiencies_BW'].loc[effi, run])})

        for procapa in list(Controller.Parameters['NetProcessEfficiencies'].index):
            if procapa not in self.ProcessCapacities.Elements:
                self.ProcessCapacities.Elements.update({procapa: {}})
            self.ProcessCapacities.Elements[procapa].update({'ID': procapa})
            for run in list(Controller.Parameters['NetProcessEfficiencies']):
                self.ProcessCapacities.Elements[procapa].update(
                    {run: json.dumps(Controller.Parameters['NetProcessEfficiencies'].loc[procapa, run])})

        for compcap in list(Controller.Parameters['CompartmentCapacities'].index):
            if compcap not in self.CompartmentCapacities.Elements:
                self.CompartmentCapacities.Elements.update({compcap: {}})
            self.CompartmentCapacities.Elements[compcap].update({'ID': compcap})
            for run in list(Controller.Parameters['CompartmentCapacities']):
                self.CompartmentCapacities.Elements[compcap].update(
                    {run: json.dumps(Controller.Parameters['CompartmentCapacities'].loc[compcap, run])})

    def exportSBtab(self, filename_SBtab):

        EnzymeCapacitiesTable_FW = self.EnzymeCapacities_FW.toSBtab(
            table_id='enzyme_forward_capacity', table_type='QuantityMatrix', table_name='Enzyme forward capacity')
        EnzymeCapacitiesTable_BW = self.EnzymeCapacities_BW.toSBtab(
            table_id='enzyme_backward_capacity', table_type='QuantityMatrix', table_name='Enzyme backward capacity')
        ProcessCapacitiesTable = self.ProcessCapacities.toSBtab(
            table_id='Machine_capacity', table_type='QuantityMatrix', table_name='Machine capacity')
        CompartmentCapacitiesTable = self.CompartmentCapacities.toSBtab(
            table_id='compartment_capacity', table_type='QuantityMatrix', table_name='Compartment capacity')

        EnzymeCapacitiesTable_FW.filename = 'EnzymeForwardCapacities.tsv'
        EnzymeCapacitiesTable_FW.change_attribute('QuantityType', 'enzyme_capacity')
        EnzymeCapacitiesTable_FW.change_attribute('Unit', '1/h')
        EnzymeCapacitiesTable_FW.change_attribute(
            'Text', 'Enzyme forward capacities (relating enzyme concentrations to reaction rates) may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        EnzymeCapacitiesTable_FW.change_attribute('TableID', 'enzyme_forward_capacity')
        EnzymeCapacitiesTable_FW.unset_attribute('Date')
        EnzymeCapacitiesTable_FW.unset_attribute('SBtabVersion')

        EnzymeCapacitiesTable_BW.filename = 'EnzymeBackwardCapacities.tsv'
        EnzymeCapacitiesTable_BW.change_attribute('QuantityType', 'enzyme_capacity')
        EnzymeCapacitiesTable_BW.change_attribute('Unit', '1/h')
        EnzymeCapacitiesTable_BW.change_attribute(
            'Text', 'Enzyme backward capacities (relating enzyme concentrations to reaction rates) may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        EnzymeCapacitiesTable_BW.change_attribute('TableID', 'enzyme_backward_capacity')
        EnzymeCapacitiesTable_BW.unset_attribute('Date')
        EnzymeCapacitiesTable_BW.unset_attribute('SBtabVersion')

        ProcessCapacitiesTable.filename = 'MachineCapacities.tsv'
        ProcessCapacitiesTable.change_attribute('QuantityType', 'machinery_capacity')
        ProcessCapacitiesTable.change_attribute('Unit', '1/h')
        ProcessCapacitiesTable.change_attribute(
            'Text', 'Machine capacities (relating machine concentrations to process rates)  may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        ProcessCapacitiesTable.change_attribute('TableID', 'machinery_capacity')
        ProcessCapacitiesTable.unset_attribute('Date')
        ProcessCapacitiesTable.unset_attribute('SBtabVersion')

        CompartmentCapacitiesTable.filename = 'CompartmentCapacities.tsv'
        CompartmentCapacitiesTable.change_attribute('QuantityType', 'compartment_capacity')
        CompartmentCapacitiesTable.change_attribute('Unit', 'AA-residues')
        CompartmentCapacitiesTable.change_attribute(
            'Text', 'Compartment capacities (defining the maximal macromolecule density in a compartment)  may depend on model parameters such as growth rate and will therefore vary between simulation runs.')
        CompartmentCapacitiesTable.change_attribute('TableID', 'compartment_capacity')
        CompartmentCapacitiesTable.unset_attribute('Date')
        CompartmentCapacitiesTable.unset_attribute('SBtabVersion')

        Out = SBtab.SBtabDocument(name='RBAparameters', sbtab_init=None,
                                  filename=str(filename_SBtab+'.tsv'))

        Out.add_sbtab(EnzymeCapacitiesTable_FW)
        Out.add_sbtab(EnzymeCapacitiesTable_BW)
        Out.add_sbtab(ProcessCapacitiesTable)
        Out.add_sbtab(CompartmentCapacitiesTable)

        Out.change_attribute('DocumentName', 'RBA parameters')
        Out.name = filename_SBtab
        Out.change_attribute('DocumentType', 'rba-simulation-parameters')
        Out.write()
