# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import sys
import os.path
import numpy
import pandas
import copy
import json
import jxmlease

# import xml.etree.ElementTree as ET
import libsbml
# package imports
import rba
from rba.core.constraint_blocks import ConstraintBlocks

from rbatools.infoMatrices import InfoMatrices
from rbatools.description_block import DescriptionBlock
from rbatools.metabolite_block import MetaboliteBlock
from rbatools.module_block import ModuleBlock
from rbatools.process_block import ProcessBlock
from rbatools.reaction_block import ReactionBlock
from rbatools.enzyme_block import EnzymeBlock
from rbatools.protein_block import ProteinBlock
from rbatools.macromolecule_block import MacromoleculeBlock
from rbatools.compartment_block import CompartmentBlock
from rbatools.metabolite_constraints import MetaboliteConstraints
from rbatools.density_constraints import DensityConstraints
from rbatools.process_constraints import ProcessConstraints
from rbatools.enzyme_constraints import EnzymeConstraints
from rbatools.statistics_block import StatisticsBlock
from rbatools.target_block import TargetBlock

from sbtab import SBtab


class RBA_ModelStructure(object):
    """
    Class holding information on model-structure.

    Attributes
    ----------
    GeneralInfo : rbatools.description_block.Description_block
         Model description
    MetaboliteInfo : rbatools.metabolite_block.Metabolite_block
         Metabolite information
    ModuleInfo : rbatools.module_block.Module_block
         Module information
    ProcessInfo : rbatools.process_block.Process_block
         Process information
    ReactionInfo : rbatools.reaction_block.Reaction_block
         Reaction information
    EnzymeInfo : rbatools.enzyme_block.Enzyme_block
         Enzyme information
    ProteinInfo : rbatools.protein_block.Protein_block
         Protein information
    MacromoleculeInfo : rbatools.macromolecule_block.MacromoleculeBlock
         Macromolecule (other than protein) information
    CompartmentInfo : rbatools.compartment_block.Compartment_block
         Compartment information
    MetaboliteConstraintsInfo : rbatools.metabolite_constraints.Metabolite_constraints
         Metabolite-constraint information
    DensityConstraintsInfo : rbatools.density_constraints.Density_constraints
         Density-constraint information
    ProcessConstraintsInfo : rbatools.process_constraints.Process_constraints
         Process-constraint information
    EnzymeConstraintsInfo : rbatools.enzyme_constraints.Enzyme_constraints
         Enzyme-constraint information
    ModelStatistics : rbatools.statistics_block.Statistics_block
         Statistics on Model
    TargetInfo : rbatools.target_block.TargetBlock
         Target information
    ProteinMatrix : numpy.array
        Matrix mapping proteins to consumers (enzymes and process-machineries)
    MediumDependencies : dict
        Dictionary with boundary metabolites as keys
        and a list of constraint_IDs, which are affected by its concentration, as values.
    MuDependencies : list
        List of constraint_IDs, which are affected by the growth-rate

    Methods
    ----------
    fromFiles:
        Generates model-structure object from model-files and provided auxilliary information.
    fromJSON:
        Generates model-structure object from model-structure in JSON-format.
    exportJSON:
        Saves model-structure object in JSON-format in specified directory, under the name ModelStructure.json.
    exportSBtab:
        Saves model-structure object in SBtab-format under the specified name (filename).
    generateMatrices:
        Generates information-matricx object and stores it as attribute 'InfoMatrices'
    """

    def fromFiles(self, xml_dir):
        """
        Generates model-structure object from model-files and provided auxilliary information.

        Parameters
        ----------
        xml_dir : str
            Directory, where RBA-model is located
        """

        UniprotFile = importUniprotFile(xml_dir)
        GeneMap = importGeneAnnotations(xml_dir)
        Info = importModelInfo(xml_dir)
        SBMLfile = str('Not There')
        if Info['Value']['SBML-file'] != 'Not Provided':
            SBMLfile = importSbmlFile(xml_dir, str(Info['Value']['SBML-file']))

        MetaboliteAnnotations = importMetaboliteAnnotations(xml_dir)
        ReactionAnnotations = importReactionAnnotations(xml_dir)

        print('')
        print('Generating model-structure')
        print('...')

        model = rba.RbaModel.from_xml(xml_dir)
        Zero_matrix = rba.ConstraintMatrix(model)
        Zero_matrix.build_matrices(0)
        constraints = sortConstraints(Zero_matrix, model)

        MetaboliteInfo = MetaboliteBlock()
        ModuleInfo = ModuleBlock()
        ProcessInfo = ProcessBlock()
        ReactionInfo = ReactionBlock()
        EnzymeInfo = EnzymeBlock()
        ProteinInfo = ProteinBlock()
        MacromoleculeInfo = MacromoleculeBlock()
        CompartmentInfo = CompartmentBlock()
        TargetInfo = TargetBlock()

        TargetInfo.fromFiles(model)
        MetaboliteInfo.fromFiles(model, Info, MetaboliteAnnotations, SBMLfile)
        ModuleInfo.fromFiles(model, SBMLfile)
        ProcessInfo.fromFiles(model, Info)
        ReactionInfo.fromFiles(model, Info, ReactionAnnotations, SBMLfile, MetaboliteInfo)
        EnzymeInfo.fromFiles(model, Info)
        ProteinInfo.fromFiles(model, GeneMap, Info, UniprotFile)
        MacromoleculeInfo.fromFiles(model)
        CompartmentInfo.fromFiles(model, Info)

        self.GeneralInfo = DescriptionBlock()
        self.GeneralInfo.fromFiles(Info)
        self.MetaboliteInfo = MetaboliteInfo
        self.ModuleInfo = ModuleInfo
        self.ProcessInfo = ProcessInfo
        self.TargetInfo = TargetInfo

        for protein in ProteinInfo.Elements.keys():
            AssociatedEnzyme = findAssociatedEnzyme(EnzymeInfo.Elements, protein)
            ProteinInfo.Elements[protein]['associatedEnzymes'] = AssociatedEnzyme['Enz']
            ProteinInfo.Elements[protein]['associatedReactions'] = AssociatedEnzyme['Rx']

        CB = ConstraintBlocks(model)
        for enzyme in EnzymeInfo.Elements.keys():
            EnzymeInfo.Elements[enzyme]['Isozymes'] = findIsozymes(
                enzyme, CB, ReactionInfo.Elements, EnzymeInfo.Elements[enzyme]['Reaction'])

        for rx in ReactionInfo.Elements.keys():
            if ReactionInfo.Elements[rx]['Enzyme'] is not '':
                ReactionInfo.Elements[rx]['Compartment_Machinery'] = EnzymeInfo.Elements[ReactionInfo.Elements[rx]
                                                                                         ['Enzyme']]['EnzymeCompartment']

        for comp in CompartmentInfo.Elements.keys():
            ContainedMacromolecules = findContainedMacromolecules(comp, MacromoleculeInfo.Elements)
            ContainedProteins = findContainedProteins(comp, ProteinInfo.Elements)
            ContainedEnzymes = findContainedEnzymes(
                ContainedProteins, ProteinInfo.Elements, EnzymeInfo.Elements)
            CompartmentInfo.Elements[comp]['associatedMacromolecules'] = ContainedMacromolecules
            CompartmentInfo.Elements[comp]['associatedProteins'] = ContainedProteins
            CompartmentInfo.Elements[comp]['associatedReactions'] = ContainedEnzymes['containedReactions']
            CompartmentInfo.Elements[comp]['associatedEnzymes'] = ContainedEnzymes['containedEnzymes']

        self.ReactionInfo = ReactionInfo
        self.EnzymeInfo = EnzymeInfo
        self.ProteinInfo = ProteinInfo
        self.MacromoleculeInfo = MacromoleculeInfo
        self.CompartmentInfo = CompartmentInfo

        self.MetaboliteConstraintsInfo = MetaboliteConstraints()
        self.DensityConstraintsInfo = DensityConstraints()
        self.ProcessConstraintsInfo = ProcessConstraints()
        self.EnzymeConstraintsInfo = EnzymeConstraints()

        self.MetaboliteConstraintsInfo.fromFiles(constraints, Zero_matrix)
        self.DensityConstraintsInfo.fromFiles(model, constraints, Zero_matrix)
        self.ProcessConstraintsInfo.fromFiles(model, constraints, Zero_matrix)
        self.EnzymeConstraintsInfo.fromFiles(model, constraints, Zero_matrix)

        BioConstraintStats = StatsConstraintsBiological(self.MetaboliteConstraintsInfo.Elements,
                                                        self.EnzymeConstraintsInfo.Elements,
                                                        self.DensityConstraintsInfo.Elements,
                                                        self.ProcessConstraintsInfo.Elements)

        MathConstraintStats = StatsConstraintsMathematical(self.MetaboliteConstraintsInfo.Elements,
                                                           self.EnzymeConstraintsInfo.Elements,
                                                           self.DensityConstraintsInfo.Elements,
                                                           self.ProcessConstraintsInfo.Elements,
                                                           self.EnzymeInfo.Elements,
                                                           self.ReactionInfo.Elements,
                                                           self.ProcessInfo.Elements)

        FullOverview = generateOverview(self.ReactionInfo.overview(),
                                        self.MetaboliteInfo.overview(),
                                        self.ModuleInfo.overview(),
                                        self.EnzymeInfo.overview(),
                                        self.ProteinInfo.overview(),
                                        self.MacromoleculeInfo.overview(),
                                        self.ProcessInfo.overview(),
                                        self.CompartmentInfo.overview(),
                                        self.TargetInfo.overview(),
                                        BioConstraintStats,
                                        MathConstraintStats)

        self.ModelStatistics = StatisticsBlock()
        self.ModelStatistics.derive(FullOverview)

        self.ProteinMatrix = generateProteinMatrix(self)
        self.ProteinGeneMatrix = generateProteinGeneMatrix(self)
        self.MediumDependencies, self.MuDependencies = findParameterDependencies(self)

    def fromJSON(self, inputString):
        """
        Generates model-structure object from model-structure in JSON-format.

        Parameters
        ----------
        inputString : json-str
            JSON-string to be parsed in to ModelStructure-object.
        """
        Block = json.loads(inputString)
        self.ModelStatistics = StatisticsBlock()
        self.GeneralInfo = DescriptionBlock()
        self.ProcessInfo = ProcessBlock()
        self.CompartmentInfo = CompartmentBlock()
        self.MetaboliteInfo = MetaboliteBlock()
        self.TargetInfo = TargetBlock()
        self.ModuleInfo = ModuleBlock()
        self.EnzymeInfo = EnzymeBlock()
        self.ProteinInfo = ProteinBlock()
        self.MacromoleculeInfo = MacromoleculeBlock()
        self.ReactionInfo = ReactionBlock()
        self.DensityConstraintsInfo = DensityConstraints()
        self.ProcessConstraintsInfo = ProcessConstraints()
        self.MetaboliteConstraintsInfo = MetaboliteConstraints()
        self.EnzymeConstraintsInfo = EnzymeConstraints()

        self.ModelStatistics.fromDict(Block['ModelStatistics'])
        self.GeneralInfo.fromDict(Block['ModelInformation'])
        self.ProcessInfo.fromDict(Block['Process'])
        self.CompartmentInfo.fromDict(Block['Compartment'])
        self.MetaboliteInfo.fromDict(Block['Metabolite'])
        self.ModuleInfo.fromDict(Block['Module'])
        self.EnzymeInfo.fromDict(Block['Enzyme'])
        self.ProteinInfo.fromDict(Block['Protein'])
        self.MacromoleculeInfo.fromDict(Block['Macromolecule'])
        self.ReactionInfo.fromDict(Block['Reaction'])
        self.DensityConstraintsInfo.fromDict(Block['DensityConstraint'])
        self.ProcessConstraintsInfo.fromDict(Block['ProcessConstraint'])
        self.MetaboliteConstraintsInfo.fromDict(Block['MetaboliteConstraint'])
        self.EnzymeConstraintsInfo.fromDict(Block['EnzymeConstraint'])
        self.TargetInfo.fromDict(Block['Target'])
        self.ProteinMatrix = Block['ProteinMatrix']
        self.ProteinMatrix['Matrix'] = numpy.array(self.ProteinMatrix['Matrix'])
        self.ProteinGeneMatrix = Block['ProteinGeneMatrix']
        self.ProteinGeneMatrix['Matrix'] = numpy.array(self.ProteinGeneMatrix['Matrix'])
        self.MediumDependencies = Block['MediumDependencies']
        self.MuDependencies = Block['MuDependencies']

    def exportJSON(self, path):
        """
        Saves model-structure object in JSON-format in specified directory, under the name ModelStructure.json.

        Parameters
        ----------
        path : str
            Directory, where to save JSON-file
        """
        Block = {'ModelInformation': self.GeneralInfo.Elements,
                 'ModelStatistics': self.ModelStatistics.Elements,
                 'Process': self.ProcessInfo.Elements,
                 'Compartment': self.CompartmentInfo.Elements,
                 'Metabolite': self.MetaboliteInfo.Elements,
                 'Target': self.TargetInfo.Elements,
                 'Module': self.ModuleInfo.Elements,
                 'Enzyme': self.EnzymeInfo.Elements,
                 'Protein': self.ProteinInfo.Elements,
                 'Macromolecule': self.MacromoleculeInfo.Elements,
                 'Reaction': self.ReactionInfo.Elements,
                 'DensityConstraint': self.DensityConstraintsInfo.Elements,
                 'ProcessConstraint': self.ProcessConstraintsInfo.Elements,
                 'MetaboliteConstraint': self.MetaboliteConstraintsInfo.Elements,
                 'EnzymeConstraint': self.EnzymeConstraintsInfo.Elements,
                 'ProteinMatrix': self.ProteinMatrix,
                 'ProteinGeneMatrix': self.ProteinGeneMatrix,
                 'MediumDependencies': self.MediumDependencies,
                 'MuDependencies': self.MuDependencies}
        # blocks.metabolism.external, model.metabolism.reactions
        # Block['ExternalStoichiometryMatrix']['S'] = Block['ProteinMatrix']['Matrix'].tolist()
        # Block['ExternalStoichiometryMatrix']['Reactions'] = Block['ProteinMatrix']['Matrix'].tolist()
        # Block['ExternalStoichiometryMatrix']['Metabolites'] = [i.id for i in blocks.metabolism.external]
        Block['ProteinMatrix']['Matrix'] = Block['ProteinMatrix']['Matrix'].tolist()
        Block['ProteinGeneMatrix']['Matrix'] = Block['ProteinGeneMatrix']['Matrix'].tolist()
        JSONstring = json.dumps(Block, default=JSON_Int64_compensation)
        filename = path + '/ModelStructure.json'
        f = open(filename, 'w')
        f.write(JSONstring)
        f.close()
        return(JSONstring)

    def exportSBtab(self, add_links=False, filename=None):
        """
        Saves model-structure object in SBtab-format under the specified name (filename).

        Parameters
        ----------
        filename : str
            Name, under which to save SBtab-file
        add_links : str
            Wheter to implement entry-format, which allows links between table-elements.
        """
        GeneralInfoTable = self.GeneralInfo.toSBtab(table_id='ModelMetadata', table_type='Quantity', table_name='Model Information', Col_list=[
                                                    'Measure', 'Value'], NameList=['Element', 'Value'])
        GeneralInfoTable.filename = 'ModelMetadata.tsv'
        GeneralInfoTable.change_attribute('Text', 'Model metadata')
        GeneralInfoTable.unset_attribute('Date')
        GeneralInfoTable.unset_attribute('SBtabVersion')

        StatsTable = self.ModelStatistics.toSBtab(
            table_id='ModelElements', table_type='Quantity', table_name='Model Size', Col_list=['Measure', 'Value'], NameList=['Element', 'Number'])
        StatsTable.filename = 'ModelElements.tsv'
        StatsTable.change_attribute('Text', 'Numbers of cell elements and constraints')
        StatsTable.unset_attribute('Date')
        StatsTable.unset_attribute('SBtabVersion')

        MetaboliteBlock_forChanges = copy.deepcopy(self.MetaboliteInfo)
        if add_links:
            for k in list(MetaboliteBlock_forChanges.Elements.keys()):
                if MetaboliteBlock_forChanges.Elements[k]['Compartment'] in self.CompartmentInfo.Elements.keys():
                    MetaboliteBlock_forChanges.Elements[k]['Compartment'] = '(!' + 'Compartment' + \
                        '/'+MetaboliteBlock_forChanges.Elements[k]['Compartment']+'!)'
                for index, id in enumerate(MetaboliteBlock_forChanges.Elements[k]['ReactionsInvolvedWith']):
                    MetaboliteBlock_forChanges.Elements[k]['ReactionsInvolvedWith'][
                        index] = '(!' + 'Reaction'+'/'+id+'!)'
                for id in list(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    if not pandas.isna(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id]):
                        MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id] = str('(!identifiers:'+id+'/'+str(
                            MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id])+'|'+str(MetaboliteBlock_forChanges.Elements[k]['OtherIDs'][id])+'!)')

        MetaboliteTable = MetaboliteBlock_forChanges.toSBtab(table_id='Compound', table_type='Quantity', table_name='Metabolites', Col_list=[
            'ID', 'Name', 'Compartment', 'Type', 'ReactionsInvolvedWith', 'OtherIDs'], NameList=['ID', 'Name', 'Compartment', 'Type', 'Reactions', 'Annotation'])
        MetaboliteTable.filename = 'Compound.tsv'
        MetaboliteTable.change_attribute(
            'Text', 'Metabolite species are localised in cell compartments and are associated with metabolite mass balance constraints.')
        MetaboliteTable.unset_attribute('Date')
        MetaboliteTable.unset_attribute('SBtabVersion')

        ReactionBlock_forChanges = copy.deepcopy(self.ReactionInfo)
        if add_links:
            for k in list(ReactionBlock_forChanges.Elements.keys()):
                oldEnz = ReactionBlock_forChanges.Elements[k]['Enzyme']
                if len(oldEnz) > 0:
                    ReactionBlock_forChanges.Elements[k]['Enzyme'] = '(!' + \
                        'Enzyme'+'/'+oldEnz+'!)'

                replacements = {}
                for i in ReactionBlock_forChanges.Elements[k]['Formula'].split(' '):
                    if i in self.MetaboliteInfo.Elements.keys():
                        replacements[i] = '(!'+'Compound' + '/'+i+'!)'
                for i in replacements.keys():
                    ReactionBlock_forChanges.Elements[k]['Formula'] = ReactionBlock_forChanges.Elements[k]['Formula'].replace(
                        str(' '+i+' '), str(' '+replacements[i]+' '))

                for index, id in enumerate(ReactionBlock_forChanges.Elements[k]['Compartment_Species']):
                    if id in self.CompartmentInfo.Elements.keys():
                        ReactionBlock_forChanges.Elements[k]['Compartment_Species'][
                            index] = '(!' + 'Compartment'+'/'+id+'!)'

                for id in list(ReactionBlock_forChanges.Elements[k]['Reactants'].keys()):
                    ReactionBlock_forChanges.Elements[k]['Reactants']['(!'+'Compound' +
                                                                      '/'+id+'!)'] = ReactionBlock_forChanges.Elements[k]['Reactants'].pop(id)
                for id in list(ReactionBlock_forChanges.Elements[k]['Products'].keys()):
                    ReactionBlock_forChanges.Elements[k]['Products']['(!'+'Compound' +
                                                                     '/'+id+'!)'] = ReactionBlock_forChanges.Elements[k]['Products'].pop(id)
                for index, id in enumerate(ReactionBlock_forChanges.Elements[k]['Twins']):
                    ReactionBlock_forChanges.Elements[k]['Twins'][index] = '(!' + \
                        'Reaction'+'/'+id+'!)'
                for index, id in enumerate(ReactionBlock_forChanges.Elements[k]['Compartment_Machinery']):
                    ReactionBlock_forChanges.Elements[k]['Compartment_Machinery'][
                        index] = '(!' + 'Compartment'+'/'+id+'!)'
                if 'ProtoID' in list(ReactionBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    ReactionBlock_forChanges.Elements[k]['OtherIDs'].pop('ProtoID')
                for id in list(ReactionBlock_forChanges.Elements[k]['OtherIDs'].keys()):
                    if not pandas.isna(ReactionBlock_forChanges.Elements[k]['OtherIDs'][id]):
                        ReactionBlock_forChanges.Elements[k]['OtherIDs'][id] = str('(!identifiers:'+id+'/'+str(
                            ReactionBlock_forChanges.Elements[k]['OtherIDs'][id])+'|'+str(ReactionBlock_forChanges.Elements[k]['OtherIDs'][id])+'!)')

        ReactionTable = ReactionBlock_forChanges.toSBtab(table_id='Reaction', table_type='Quantity', table_name='Reactions', Col_list=['ID', 'Name', 'Type', 'Reversible', 'Formula', 'Enzyme', 'Compartment_Machinery', 'Twins', 'Compartment_Species', 'OtherIDs'], NameList=[
            'ID', 'Name', 'Type', 'IsReversible', 'ReactionFormula', 'Enzyme', 'EnzymeCompartment', 'IsoenzymeReactions', 'ReactionCompartment', 'Annotation'])
        ReactionTable.filename = 'Reaction.tsv'
        ReactionTable.change_attribute(
            'Text', 'Chemical reactions are localised in cell compartments. All reactions are enzyme catalysed.')
        ReactionTable.unset_attribute('Date')
        ReactionTable.unset_attribute('SBtabVersion')

        EnzymeBlock_forChanges = copy.deepcopy(self.EnzymeInfo)
        if add_links:
            for k in list(EnzymeBlock_forChanges.Elements.keys()):
                oldRx = EnzymeBlock_forChanges.Elements[k]['Reaction']
                EnzymeBlock_forChanges.Elements[k]['Reaction'] = '(!'+'Reaction'+'/'+oldRx+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['Isozymes']):
                    EnzymeBlock_forChanges.Elements[k]['Isozymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['IdenticalEnzymes']):
                    EnzymeBlock_forChanges.Elements[k]['IdenticalEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(EnzymeBlock_forChanges.Elements[k]['EnzymeCompartment']):
                    EnzymeBlock_forChanges.Elements[k]['EnzymeCompartment'][index] = '(!' + \
                        'Compartment'+'/' + id+'!)'
                for id in list(EnzymeBlock_forChanges.Elements[k]['Subunits'].keys()):
                    EnzymeBlock_forChanges.Elements[k]['Subunits']['(!'+'Protein'+'/' +
                                                                   id+'!)'] = EnzymeBlock_forChanges.Elements[k]['Subunits'].pop(id)

        EnzymeTable = EnzymeBlock_forChanges.toSBtab(table_id='Enzyme', table_type='Quantity', table_name='Enzymes', Col_list=[
            'ID', 'Reaction', 'IdenticalEnzymes', 'Subunits', 'EnzymeCompartment', 'OtherIDs', 'Isozymes'], NameList=['ID', 'CatalysedReaction', 'OtherEnzymaticActivities', 'Subunits', 'Compartment', 'Annotation', 'Isoenzymes'])
        EnzymeTable.filename = 'Enzyme.tsv'
        EnzymeTable.change_attribute(
            'Text', 'Enzymes are localised in cell compartments and catalyse specific reactions. To describe multi-functional enzyme, RBA uses multiple enzyme entries. Enzymes are associated with enzyme capacity constraints.')
        EnzymeTable.unset_attribute('Date')
        EnzymeTable.unset_attribute('SBtabVersion')

        ProteinBlock_forChanges = copy.deepcopy(self.ProteinInfo)
        if add_links:
            for k in list(ProteinBlock_forChanges.Elements.keys()):
                oldComp = ProteinBlock_forChanges.Elements[k]['Compartment']
                ProteinBlock_forChanges.Elements[k]['Compartment'] = '(!' + \
                    'Compartment'+'/'+oldComp+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['associatedReactions']):
                    ProteinBlock_forChanges.Elements[k]['associatedReactions'][
                        index] = '(!' + 'Reaction'+'/'+id+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['associatedEnzymes']):
                    ProteinBlock_forChanges.Elements[k]['associatedEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(ProteinBlock_forChanges.Elements[k]['SupportsProcess']):
                    ProteinBlock_forChanges.Elements[k]['SupportsProcess'][index] = '(!' + \
                        'Process'+'/'+id+'!)'
                for id in list(ProteinBlock_forChanges.Elements[k]['ProcessRequirements'].keys()):
                    ProteinBlock_forChanges.Elements[k]['ProcessRequirements'][
                        '(!'+'Process'+'/' + id+'!)'] = ProteinBlock_forChanges.Elements[k]['ProcessRequirements'].pop(id)
                for id in list(ProteinBlock_forChanges.Elements[k]['ExternalIDs'].keys()):
                    if not pandas.isna(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id]):
                        if id == 'UniprotID':
                            ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] = str('(!identifiers:uniprot/'+str(
                                ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'|'+str(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'!)')
                        if id == 'ECnumber':
                            if ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] != 'nan':
                                ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id] = str('(!identifiers:ec-code/'+str(
                                    ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id]).split('EC ')[1]+'|'+str(ProteinBlock_forChanges.Elements[k]['ExternalIDs'][id])+'!)')

        ProteinTable = ProteinBlock_forChanges.toSBtab(table_id='Protein', table_type='Quantity', table_name='Proteins', Col_list=['ID', 'ProtoID', 'Name', 'Compartment', 'SupportsProcess', 'associatedEnzymes', 'associatedReactions', 'AAnumber', 'ProcessRequirements', 'Weight', 'ExternalIDs', 'Function'], NameList=[
            'ID', 'ProtoID', 'Name', 'Compartment', 'ContributesToProcess', 'EnzymaticActivity', 'CatalysedReaction', 'ChainLength', 'RequiresProcess', 'MolecularWeight', 'Annotation', 'Function'])
        ProteinTable.filename = 'Protein.tsv'
        ProteinTable.change_attribute(
            'Text', 'Proteins are localised in cell compartments and can be (or be part of) metabolic enzymes.')
        ProteinTable.unset_attribute('Date')
        ProteinTable.unset_attribute('SBtabVersion')

        MacromoleculeBlock_forChanges = copy.deepcopy(self.MacromoleculeInfo)
        if add_links:
            for k in list(MacromoleculeBlock_forChanges.Elements.keys()):
                oldComp = MacromoleculeBlock_forChanges.Elements[k]['Compartment']
                MacromoleculeBlock_forChanges.Elements[k]['Compartment'] = '(!' + \
                    'Compartment'+'/'+oldComp+'!)'
                for index, id in enumerate(MacromoleculeBlock_forChanges.Elements[k]['SupportsProcess']):
                    MacromoleculeBlock_forChanges.Elements[k][
                        'SupportsProcess'][index] = '(!' + 'Process'+'/'+id+'!)'
                for id in list(MacromoleculeBlock_forChanges.Elements[k]['ProcessRequirements'].keys()):
                    MacromoleculeBlock_forChanges.Elements[k]['ProcessRequirements'][
                        '(!'+'Process'+'/' + id+'!)'] = MacromoleculeBlock_forChanges.Elements[k]['ProcessRequirements'].pop(id)

        MacromoleculeTable = MacromoleculeBlock_forChanges.toSBtab(table_id='Macromolecule', table_type='Quantity', table_name='Macromolecules', Col_list=[
                                                                   'ID', 'ProtoID', 'Type', 'Compartment', 'SupportsProcess', 'ProcessRequirements'], NameList=['ID', 'ProtoID', 'Type', 'Compartment', 'SupportedtedProcess', 'RequiresProcess'])
        MacromoleculeTable.filename = 'Macromolecules.tsv'
        MacromoleculeTable.change_attribute(
            'Text', 'Macromolecule are localised in cell compartments and can be part of cellular machinery or serve another function inside the cell.')
        MacromoleculeTable.unset_attribute('Date')
        MacromoleculeTable.unset_attribute('SBtabVersion')

        CompartmentBlock_forChanges = copy.deepcopy(self.CompartmentInfo)
        if add_links:
            for k in list(CompartmentBlock_forChanges.Elements.keys()):
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedReactions']):
                    CompartmentBlock_forChanges.Elements[k]['associatedReactions'][
                        index] = '(!' + 'Reaction'+'/'+id+'!)'
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedEnzymes']):
                    CompartmentBlock_forChanges.Elements[k]['associatedEnzymes'][index] = '(!' + \
                        'Enzyme'+'/'+id+'!)'
                for index, id in enumerate(CompartmentBlock_forChanges.Elements[k]['associatedProteins']):
                    CompartmentBlock_forChanges.Elements[k]['associatedProteins'][
                        index] = '(!' + 'Protein'+'/'+id+'!)'

        CompartmentTable = CompartmentBlock_forChanges.toSBtab(table_id='Compartment', table_type='Quantity', table_name='Cell Compartments', Col_list=[
                                                               'ID', 'associatedProteins', 'associatedEnzymes', 'associatedReactions'], NameList=['ID', 'Proteins', 'Enzymes', 'Reactions'])

        CompartmentTable.filename = 'Compartment.tsv'
        CompartmentTable.change_attribute(
            'Text', 'Cell compartments are used to describe the localisation of proteins, enzymes, and reactions and are associated with density constraints.')
        CompartmentTable.unset_attribute('Date')
        CompartmentTable.unset_attribute('SBtabVersion')

        ModuleBlock_forChanges = copy.deepcopy(self.ModuleInfo)
        if add_links:
            for k in list(ModuleBlock_forChanges.Elements.keys()):
                for index, id in enumerate(ModuleBlock_forChanges.Elements[k]['Members']):
                    if id in self.ReactionInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Reaction'+'/'+id+'!)'
                    if id in self.ProteinInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Protein'+'/'+id+'!)'
                    if id in self.MetaboliteInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Compound'+'/'+id+'!)'
                    if id in self.CompartmentInfo.Elements.keys():
                        ModuleBlock_forChanges.Elements[k]['Members'][index] = '(!' + \
                            'Compartment'+'/'+id+'!)'

        ModuleTable = ModuleBlock_forChanges.toSBtab(table_id='CellModule', table_type='Quantity', table_name='Cell Modules', Col_list=[
                                                     'ID', 'Name', 'Members'], NameList=['ID', 'Name', 'Contains'])
        ModuleTable.filename = 'Module.tsv'
        ModuleTable.change_attribute('Text', 'Information on Modules in RBAmodel')
        ModuleTable.unset_attribute('Date')
        ModuleTable.unset_attribute('SBtabVersion')

        ProcessBlock_forChanges = copy.deepcopy(self.ProcessInfo)
        if add_links:
            for k in list(ProcessBlock_forChanges.Elements.keys()):
                for component in ProcessBlock_forChanges.Elements[k]['Components'].keys():
                    ProductKeys = list(
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Products'].keys())
                    for species in ProductKeys:
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Products'][
                            '(!'+'Compound'+'/'+species+'!)'] = ProcessBlock_forChanges.Elements[k]['Components'][component]['Products'].pop(species)
                    ReactantKeys = list(
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Reactants'].keys())
                    for species in ReactantKeys:
                        ProcessBlock_forChanges.Elements[k]['Components'][component]['Reactants'][
                            '(!'+'Compound'+'/'+species+'!)'] = ProcessBlock_forChanges.Elements[k]['Components'][component]['Reactants'].pop(species)
                for id in list(ProcessBlock_forChanges.Elements[k]['Composition'].keys()):
                    if id in self.ProteinInfo.Elements.keys():
                        ProcessBlock_forChanges.Elements[k]['Composition']['(!'+'Protein'+'/' +
                                                                           id+'!)'] = ProcessBlock_forChanges.Elements[k]['Composition'].pop(id)
                    elif id in self.MacromoleculeInfo.Elements.keys():
                        ProcessBlock_forChanges.Elements[k]['Composition']['(!'+'Macromolecule' +
                                                                           '/' + id+'!)'] = ProcessBlock_forChanges.Elements[k]['Composition'].pop(id)
                    elif id in self.MetaboliteInfo.Elements.keys():
                        ProcessBlock_forChanges.Elements[k]['Composition']['(!'+'Metabolite'+'/' +
                                                                           id+'!)'] = ProcessBlock_forChanges.Elements[k]['Composition'].pop(id)

        ProcessTable = ProcessBlock_forChanges.toSBtab(table_id='Process', table_type='Quantity', table_name='Macromolecular Processes', Col_list=[
            'ID', 'Name', 'Composition', 'Components', 'Initiation'], NameList=['ID', 'Name', 'MachineSubunits', 'MachineComponents', 'InitiationCofactors'])
        ProcessTable.filename = 'Process.tsv'
        ProcessTable.change_attribute(
            'Text', 'Macromolecular machines catalyse the biochemical reactions that produce, modify, and degrade macromolecules. They are catalysed by macromolecular machines and associated with process capacity constraints.')
        ProcessTable.unset_attribute('Date')
        ProcessTable.unset_attribute('SBtabVersion')

        MetaboliteConstraintsBlock_forChanges = copy.deepcopy(self.MetaboliteConstraintsInfo)
        if add_links:
            for k in list(MetaboliteConstraintsBlock_forChanges.Elements.keys()):
                oldMet = MetaboliteConstraintsBlock_forChanges.Elements[k]['AssociatedMetabolite']
                MetaboliteConstraintsBlock_forChanges.Elements[k][
                    'AssociatedMetabolite'] = '(!' + 'Compound'+'/'+oldMet+'!)'

        MetaboliteConstraintTable = MetaboliteConstraintsBlock_forChanges.toSBtab(table_id='MetaboliteConstraint', table_type='Quantity', table_name='Metabolite Mass Balance Constraints', Col_list=[
                                                                                  'ID', 'AssociatedMetabolite', 'Type'], NameList=['ID', 'Metabolite', 'Type'])
        MetaboliteConstraintTable.filename = 'MetaboliteConstraint.tsv'
        MetaboliteConstraintTable.change_attribute(
            'Text', 'Metabolite mass balance constraints ensure mass balance and stationary fluxes in metabolism.')
        MetaboliteConstraintTable.unset_attribute('Date')
        MetaboliteConstraintTable.unset_attribute('SBtabVersion')

        CompartmentConstraintsBlock_forChanges = copy.deepcopy(self.DensityConstraintsInfo)
        if add_links:
            for k in list(CompartmentConstraintsBlock_forChanges.Elements.keys()):
                oldComp = CompartmentConstraintsBlock_forChanges.Elements[k]['AssociatedCompartment']
                CompartmentConstraintsBlock_forChanges.Elements[k][
                    'AssociatedCompartment'] = '(!' + 'Compartment'+'/'+oldComp+'!)'

        DensityConstraintTable = CompartmentConstraintsBlock_forChanges.toSBtab(table_id='DensityConstraint', table_type='Quantity', table_name='Density Constraints', Col_list=[
                                                                                'ID', 'AssociatedCompartment', 'Type', 'CapacityParameter'], NameList=['ID', 'Compartment', 'Type', 'Formula'])
        DensityConstraintTable.filename = 'DensityConstraint.tsv'
        DensityConstraintTable.change_attribute(
            'Text', 'Density constraints put an upper bound on the sum of macromolecule concentrations in a given compartment. The capacity parameter defines this bound in units corresponding to amino acids (contained in proteins), or one third of nucleotides (contained in RNA).')
        DensityConstraintTable.unset_attribute('Date')
        DensityConstraintTable.unset_attribute('SBtabVersion')

        ProcessConstraintsBlock_forChanges = copy.deepcopy(self.ProcessConstraintsInfo)
        if add_links:
            for k in list(ProcessConstraintsBlock_forChanges.Elements.keys()):
                oldComp = ProcessConstraintsBlock_forChanges.Elements[k]['AssociatedProcess']
                ProcessConstraintsBlock_forChanges.Elements[k]['AssociatedProcess'] = '(!' + \
                    'Process'+'/'+oldComp+'!)'

        ProcessConstraintTable = ProcessConstraintsBlock_forChanges.toSBtab(table_id='MachineryCapacityConstraint', table_type='Quantity', table_name='Machinery Capacity Constraints', Col_list=[
            'ID', 'AssociatedProcess', 'Type', 'CapacityParameter'], NameList=['ID', 'Process', 'Type', 'Formula'])
        ProcessConstraintTable.filename = 'MachineryConstraint.tsv'
        ProcessConstraintTable.change_attribute(
            'Text', 'A machinery capacity constraint states that the rate of a macromolecular process is proportional to the concentration of the catalysing machine. The proportionality constant  (capacity parameter) can be a constant or a function of model parameters such as the growth rate.')
        ProcessConstraintTable.unset_attribute('Date')
        ProcessConstraintTable.unset_attribute('SBtabVersion')

        EnzymeConstraintsBlock_forChanges = copy.deepcopy(self.EnzymeConstraintsInfo)
        if add_links:
            for k in list(EnzymeConstraintsBlock_forChanges.Elements.keys()):
                oldComp = EnzymeConstraintsBlock_forChanges.Elements[k]['AssociatedEnzyme']
                EnzymeConstraintsBlock_forChanges.Elements[k]['AssociatedEnzyme'] = '(!' + \
                    'Enzyme'+'/'+oldComp+'!)'

        EnzymeConstraintTable = EnzymeConstraintsBlock_forChanges.toSBtab(table_id='EnzymeCapacityConstraint', table_type='Quantity', table_name='Enzyme Capacity Constraints', Col_list=[
            'ID', 'AssociatedEnzyme', 'AssociatedReaction', 'Direction', 'Type', 'CapacityParameter'], NameList=['ID', 'Enzyme', 'Reaction', 'Direction', 'Type', 'Formula'])
        EnzymeConstraintTable.filename = 'EnzymeConstraint.tsv'
        EnzymeConstraintTable.change_attribute(
            'Text', 'An enzyme capacity constraint states that a reaction rate is proportional to the concentration of the catalysing enzyme. The proportionality constant (capacity parameter) can be a constant or a function of model parameters such as the growth rate.')
        EnzymeConstraintTable.unset_attribute('Date')
        EnzymeConstraintTable.unset_attribute('SBtabVersion')

        TargetBlock_forChanges = copy.deepcopy(self.TargetInfo)
        if add_links:
            for k in list(TargetBlock_forChanges.Elements.keys()):
                oldTargSpec = TargetBlock_forChanges.Elements[k]['TargetEntity']
                if oldTargSpec in list(self.MetaboliteInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Compound'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.ReactionInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Reaction'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.ProteinInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Protein'+'/'+oldTargSpec+'!)'
                if oldTargSpec in list(self.MacromoleculeInfo.Elements.keys()):
                    TargetBlock_forChanges.Elements[k]['TargetEntity'] = '(!' + \
                        'Macromolecule'+'/'+oldTargSpec+'!)'

        TargetTable = TargetBlock_forChanges.toSBtab(table_id='CellTarget', table_type='Quantity', table_name='Cell Targets', Col_list=[
            'ID', 'Group', 'Type', 'TargetEntity', 'TargetValue'], NameList=['ID', 'TargetGroup', 'TargetType', 'TargetEntity', 'TargetParameter'])
        TargetTable.filename = 'Target.tsv'
        TargetTable.change_attribute(
            'Text', 'Cell targets are biochemical variables (fluxes or concentrations) that need to remain below or above a certain threshold value to ensure a viable cell. They define constraints in the RBA model.')
        TargetTable.unset_attribute('Date')
        TargetTable.unset_attribute('SBtabVersion')

        if filename is not None:
            filename_SBtab = filename
        else:
            if add_links:
                filename_SBtab = 'RBA_model_withLinks'
            else:
                filename_SBtab = 'RBA_model'

        Out = SBtab.SBtabDocument(name='rbatools_withLinks', sbtab_init=None,
                                  filename=str(filename_SBtab+'.tsv'))
        Out.add_sbtab(GeneralInfoTable)
        Out.add_sbtab(StatsTable)
        Out.add_sbtab(MetaboliteTable)
        Out.add_sbtab(ReactionTable)
        Out.add_sbtab(EnzymeTable)
        Out.add_sbtab(ProteinTable)
        Out.add_sbtab(MacromoleculeTable)
        Out.add_sbtab(CompartmentTable)
        Out.add_sbtab(ModuleTable)
        Out.add_sbtab(ProcessTable)
        Out.add_sbtab(TargetTable)
        Out.add_sbtab(MetaboliteConstraintTable)
        Out.add_sbtab(DensityConstraintTable)
        Out.add_sbtab(ProcessConstraintTable)
        Out.add_sbtab(EnzymeConstraintTable)

        Out.name = filename_SBtab
        Out.change_attribute('DocumentName', self.GeneralInfo.Elements['Name'])
        Out.change_attribute('DocumentType', 'rba-model-synopsis')
        Out.write()

    def generateMatrices(self):
        """
        Generates information-matricx object and stores it as attribute 'InfoMatrices'
        """
        self.InfoMatrices = InfoMatrices(self)


def findParameterDependencies(ModelStructure):
    MedDepts = {}
    MuDepts = []
    for eK in ModelStructure.EnzymeConstraintsInfo.Elements.keys():
        e = ModelStructure.EnzymeConstraintsInfo.Elements[eK]
        for pf in e['CapacityParameter']:
            iV = list(pf.values())[0]['IndependentVariable']
            if list(pf.values())[0]['FunctionType'] != 'constant':
                if iV.startswith('M_'):
                    if iV.rsplit('_', 1)[0] in MedDepts.keys():
                        MedDepts[iV.rsplit('_', 1)[0]].append(e['ID'])
                    elif iV.rsplit('_', 1)[0] not in MedDepts.keys():
                        MedDepts.update({iV.rsplit('_', 1)[0]: [e['ID']]})
                if iV == 'growth_rate':
                    if e['ID'] not in MuDepts:
                        MuDepts.append(e['ID'])
    for dK in ModelStructure.DensityConstraintsInfo.Elements.keys():
        d = ModelStructure.DensityConstraintsInfo.Elements[dK]
        for pf in d['CapacityParameter']:
            iV = list(pf.values())[0]['IndependentVariable']
            if list(pf.values())[0]['FunctionType'] != 'constant':
                if iV.startswith('M_'):
                    if iV.rsplit('_', 1)[0] in MedDepts.keys():
                        MedDepts[iV.rsplit('_', 1)[0]].append(d['ID'])
                    elif iV.rsplit('_', 1)[0] not in MedDepts.keys():
                        MedDepts.update({iV.rsplit('_', 1)[0]: [d['ID']]})
                if iV == 'growth_rate':
                    if d['ID'] not in MuDepts:
                        MuDepts.append(d['ID'])
    for pK in ModelStructure.ProcessConstraintsInfo.Elements.keys():
        p = ModelStructure.ProcessConstraintsInfo.Elements[pK]
        for pf in p['CapacityParameter']:
            iV = list(pf.values())[0]['IndependentVariable']
            if list(pf.values())[0]['FunctionType'] != 'constant':
                if iV.startswith('M_'):
                    if iV.rsplit('_', 1)[0] in MedDepts.keys():
                        MedDepts[iV.rsplit('_', 1)[0]].append(p['ID'])
                    elif iV.rsplit('_', 1)[0] not in MedDepts.keys():
                        MedDepts.update({iV.rsplit('_', 1)[0]: [p['ID']]})
                if iV == 'growth_rate':
                    if p['ID'] not in MuDepts:
                        MuDepts.append(p['ID'])
    return(MedDepts, MuDepts)


def generateProteinGeneMatrix(ModelStructure):
    uniqueProteins = []
    uniqueProteinMap = {}
    Proteins = list(ModelStructure.ProteinInfo.Elements.keys())
    for i in ModelStructure.ProteinInfo.Elements.keys():
        if ModelStructure.ProteinInfo.Elements[i]['ProtoID'] not in list(uniqueProteinMap.keys()):
            uniqueProteinMap.update({ModelStructure.ProteinInfo.Elements[i]['ProtoID']: []})
            uniqueProteins.append(ModelStructure.ProteinInfo.Elements[i]['ProtoID'])
        uniqueProteinMap[ModelStructure.ProteinInfo.Elements[i]['ProtoID']].append(i)
    ProteinProteinMatrix = numpy.zeros(
        (len(list(uniqueProteinMap.keys())), len(list(ModelStructure.ProteinInfo.Elements.keys()))))
    for u in list(uniqueProteinMap.keys()):
        row_ind = uniqueProteins.index(u)
        for i in uniqueProteinMap[u]:
            col_ind = Proteins.index(i)
            ProteinProteinMatrix[row_ind, col_ind] = 1
    return({'Matrix': numpy.array(ProteinProteinMatrix), 'Proteins': Proteins, 'ProtoProteins': uniqueProteins})


def generateProteinMatrix(ModelStructure):
    Proteins = list(ModelStructure.ProteinInfo.Elements.keys())
    # print(list(ModelStructure.ProcessInfo.Elements.keys()))
    Processes = [ModelStructure.ProcessInfo.Elements[i]['ID'] +
                 '_machinery' for i in list(ModelStructure.ProcessInfo.Elements.keys())]
    Enzymes = list(ModelStructure.EnzymeInfo.Elements.keys())
    Consumers = list(set(list(Enzymes+Processes)))
    ProteinMatrix = numpy.zeros((len(Proteins), len(Consumers)))
    for p in Proteins:
        if len(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']) > 0:
            # print(list(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']))
            for pc in list(ModelStructure.ProteinInfo.Elements[p]['SupportsProcess']):
                coeff = 0
                row_ind = Proteins.index(p)
                col_ind = Consumers.index(
                    ModelStructure.ProcessInfo.Elements[pc]['ID']+'_machinery')
                coeff = ModelStructure.ProcessInfo.Elements[pc]['Composition'][p]
                ProteinMatrix[row_ind, col_ind] += coeff
        if len(ModelStructure.ProteinInfo.Elements[p]['associatedEnzymes']) > 0:
            for ez in list(ModelStructure.ProteinInfo.Elements[p]['associatedEnzymes']):
                coeff = 0
                row_ind = Proteins.index(p)
                col_ind = Consumers.index(ez)
                coeff = ModelStructure.EnzymeInfo.Elements[ez]['Subunits'][p]
                ProteinMatrix[row_ind, col_ind] += coeff
    return({'Matrix': numpy.array(ProteinMatrix), 'Consumers': Consumers, 'Proteins': Proteins})


def importBiggMetabolites(xml_dir):
    if os.path.isfile(str(xml_dir+'/bigg_models_metabolites.txt')):
        return(pandas.read_csv(str(xml_dir+'/bigg_models_metabolites.txt'), sep='\t', index_col=0))
    else:
        sys.exit('\n Required BiGG Metabolite File "bigg_models_metabolites.txt" not found.\n' +
                 ' Please provide in input-directory\n' +
                 ' To be found under: http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt\n')


def importBiggReactions(xml_dir):
    if os.path.isfile(str(xml_dir+'/bigg_models_reactions.txt')):
        return(pandas.read_csv(str(xml_dir+'/bigg_models_reactions.txt'), sep='\t', index_col=0))
    else:
        sys.exit('\n Required BiGG Reaction File "bigg_models_reactions.txt "not found.\n' +
                 ' Please provide in input-directory\n' +
                 ' To be found under: http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt\n')


def importUniprotFile(xml_dir):
    if os.path.isfile(str(xml_dir+'/uniprot.csv')):
        return(pandas.read_csv(str(xml_dir+'/uniprot.csv'), sep='\t'))
    else:
        print('\n Uniprot-file "uniprot.csv" not found.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importSbmlFile(xml_dir, filename):
    if os.path.isfile(str(xml_dir+'/'+filename)):
        SBfile = libsbml.readSBML(str(xml_dir+'/'+filename))
        if SBfile.getNumErrors() > 0:
            SBfile.printErrors()
            print('Invalid SBML')
            return(str('Not There'))
        else:
            sbml = SBfile
            return(sbml)
    else:
        print('\n SBML-file {} not found.\n' +
              ' Continuing without additional information...\n'.format(filename))
        return(str('Not There'))


def importGeneAnnotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/GeneAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/GeneAnnotations.csv'), sep=',', index_col=0)
        if len(list(out)) == 0:
            print('WARNING: Your file "GeneAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No Gene-annotation file "GeneAnnotations.csv" provided.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importReactionAnnotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/ReactionAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/ReactionAnnotations.csv'), sep=',', index_col=0)
        if len(list(out)) == 0:
            print('WARNING: Your file "ReactionAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No Reaction-annotation file "ReactionAnnotations.csv" provided.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importMetaboliteAnnotations(xml_dir):
    if os.path.isfile(str(xml_dir+'/MetaboliteAnnotations.csv')):
        out = pandas.read_csv(str(xml_dir+'/MetaboliteAnnotations.csv'), sep=',', index_col=0)
        if len(list(out)) == 0:
            print('WARNING: Your file "MetaboliteAnnotations.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No Reaction-annotation file "MetaboliteAnnotations.csv" provided.\n' +
              ' Continuing without additional information...\n')
        return(str('Not There'))


def importModelInfo(xml_dir):
    if os.path.isfile(str(xml_dir+'/ModelInformation.csv')):
        out = pandas.read_csv(str(xml_dir+'/ModelInformation.csv'),
                              sep=',', header=0)
        out.index = list(out['Key'])
        if len(list(out)) == 0:
            print('WARNING: Your file "ModelInformation.csv" seems to be empty or has the wrong delimiter (comma required).')
        return(out)
    else:
        print('\n No model-info file "ModelInformation.csv" provided.\n' +
              ' Using dummy-information\n')
        return(pandas.DataFrame([['Name', 'ModelName'], ['Author', 'John Doe'], ['Organism', 'Life'], ['Reconstruction', 'GSMM'], ['SBML-file', 'Not Provided']], index=['Name', 'Author', 'Organism', 'Reconstruction', 'SBML-file'], columns=['Key', 'Value']))


def htmlStyle(structOriginal):
    struct = copy.deepcopy(structOriginal)

    for j in struct.ProcessInfo.Elements.keys():
        for i in struct.ProcessInfo.Elements[j]['Composition'].keys():
            struct.ProcessInfo.Elements[j]['Composition']['Protein##' +
                                                          i] = struct.ProcessInfo.Elements[j]['Composition'].pop(i)

    for j in struct.EnzymeInfo.Elements.keys():
        for i in struct.EnzymeInfo.Elements[j]['Subunits'].keys():
            struct.EnzymeInfo.Elements[j]['Subunits']['Protein##' +
                                                      i] = struct.EnzymeInfo.Elements[j]['Subunits'].pop(i)
        for i in range(len(struct.EnzymeInfo.Elements[j]['Isozymes'])):
            struct.EnzymeInfo.Elements[j]['Isozymes'][i] = 'Enzyme##' + \
                struct.EnzymeInfo.Elements[j]['Isozymes'][i]
        for i in range(len(struct.EnzymeInfo.Elements[j]['IdenticalEnzymes'])):
            struct.EnzymeInfo.Elements[j]['IdenticalEnzymes'][i] = 'Enzyme##' + \
                struct.EnzymeInfo.Elements[j]['IdenticalEnzymes'][i]
        for i in range(len(struct.EnzymeInfo.Elements[j]['EnzymeCompartment'])):
            struct.EnzymeInfo.Elements[j]['EnzymeCompartment'][i] = 'Compartment##' + \
                struct.EnzymeInfo.Elements[j]['EnzymeCompartment'][i]
        struct.EnzymeInfo.Elements[j]['Reaction'] = 'Reaction##' + \
            struct.EnzymeInfo.Elements[j]['Reaction']

    for j in struct.ReactionInfo.Elements.keys():
        for i in struct.ReactionInfo.Elements[j]['Reactants'].keys():
            struct.ReactionInfo.Elements[j]['Reactants']['Metabolite##' +
                                                         i] = struct.ReactionInfo.Elements[j]['Reactants'].pop(i)
        for i in struct.ReactionInfo.Elements[j]['Products'].keys():
            struct.ReactionInfo.Elements[j]['Products']['Metabolite##' +
                                                        i] = struct.ReactionInfo.Elements[j]['Products'].pop(i)
        for i in range(len(struct.ReactionInfo.Elements[j]['Twins'])):
            struct.ReactionInfo.Elements[j]['Twins'][i] = 'Reaction##' + \
                struct.ReactionInfo.Elements[j]['Twins'][i]
        struct.ReactionInfo.Elements[j]['Enzyme'] = 'Enzyme##' + \
            struct.ReactionInfo.Elements[j]['Enzyme']

    for j in struct.ProteinInfo.Elements.keys():
        for i in struct.ProteinInfo.Elements[j]['ProcessRequirements'].keys():
            struct.ProteinInfo.Elements[j]['ProcessRequirements']['Process##' +
                                                                  i] = struct.ProteinInfo.Elements[j]['ProcessRequirements'].pop(i)
        for i in range(len(struct.ProteinInfo.Elements[j]['associatedReactions'])):
            struct.ProteinInfo.Elements[j]['associatedReactions'][i] = 'Reaction##' + \
                struct.ProteinInfo.Elements[j]['associatedReactions'][i]
        for i in range(len(struct.ProteinInfo.Elements[j]['associatedEnzymes'])):
            struct.ProteinInfo.Elements[j]['associatedEnzymes'][i] = 'Enzyme##' + \
                struct.ProteinInfo.Elements[j]['associatedEnzymes'][i]
        for i in range(len(struct.ProteinInfo.Elements[j]['SupportsProcess'])):
            struct.ProteinInfo.Elements[j]['SupportsProcess'][i] = 'Process##' + \
                struct.ProteinInfo.Elements[j]['SupportsProcess'][i]
        struct.ProteinInfo.Elements[j]['Compartment'] = 'Compartment##' + \
            struct.ProteinInfo.Elements[j]['Compartment']

    for j in struct.CompartmentInfo.Elements.keys():
        for i in range(len(struct.CompartmentInfo.Elements[j]['associatedProteins'])):
            struct.CompartmentInfo.Elements[j]['associatedProteins'][i] = 'Protein##' + \
                struct.CompartmentInfo.Elements[j]['associatedProteins'][i]
        for i in range(len(struct.CompartmentInfo.Elements[j]['associatedReactions'])):
            struct.CompartmentInfo.Elements[j]['associatedReactions'][i] = 'Reaction##' + \
                struct.CompartmentInfo.Elements[j]['associatedReactions'][i]
        for i in range(len(struct.CompartmentInfo.Elements[j]['associatedEnzymes'])):
            struct.CompartmentInfo.Elements[j]['associatedEnzymes'][i] = 'Enzyme##' + \
                struct.CompartmentInfo.Elements[j]['associatedEnzymes'][i]

    for j in struct.MetaboliteInfo.Elements.keys():
        for i in range(len(struct.MetaboliteInfo.Elements[j]['ReactionsInvolvedWith'])):
            struct.MetaboliteInfo.Elements[j]['ReactionsInvolvedWith'][i] = 'Reaction##' + \
                struct.MetaboliteInfo.Elements[j]['ReactionsInvolvedWith'][i]

    for j in struct.MetaboliteConstraintsInfo.Elements.keys():
        struct.MetaboliteConstraintsInfo.Elements[j]['AssociatedMetabolite'] = 'Metabolite##' + \
            struct.MetaboliteConstraintsInfo.Elements[j]['AssociatedMetabolite']
    for j in struct.EnzymeConstraintsInfo.Elements.keys():
        struct.EnzymeConstraintsInfo.Elements[j]['AssociatedEnzyme'] = 'Enzyme##' + \
            struct.EnzymeConstraintsInfo.Elements[j]['AssociatedEnzyme']
    for j in struct.EnzymeConstraintsInfo.Elements.keys():
        struct.EnzymeConstraintsInfo.Elements[j]['AssociatedReaction'] = 'Reaction##' + \
            struct.EnzymeConstraintsInfo.Elements[j]['AssociatedReaction']
    for j in struct.ProcessConstraintsInfo.Elements.keys():
        struct.ProcessConstraintsInfo.Elements[j]['AssociatedProcess'] = 'Process##' + \
            struct.ProcessConstraintsInfo.Elements[j]['AssociatedProcess']
    for j in struct.DensityConstraintsInfo.Elements.keys():
        struct.DensityConstraintsInfo.Elements[j]['AssociatedCompartment'] = 'Compartment##' + \
            struct.DensityConstraintsInfo.Elements[j]['AssociatedCompartment']

    for i in struct.CompartmentInfo.Elements.keys():
        struct.CompartmentInfo.Elements['ID_' + i] = struct.CompartmentInfo.Elements.pop(i)
    for i in struct.ProcessInfo.Elements.keys():
        struct.ProcessInfo.Elements['ID_' + i] = struct.ProcessInfo.Elements.pop(i)
    for i in struct.MetaboliteInfo.Elements.keys():
        struct.MetaboliteInfo.Elements['ID_' + i] = struct.MetaboliteInfo.Elements.pop(i)
    for i in struct.ModuleInfo.Elements.keys():
        struct.ModuleInfo.Elements['ID_' + i] = struct.ModuleInfo.Elements.pop(i)
    for i in struct.EnzymeInfo.Elements.keys():
        struct.EnzymeInfo.Elements['ID_' + i] = struct.EnzymeInfo.Elements.pop(i)
    for i in struct.ProteinInfo.Elements.keys():
        struct.ProteinInfo.Elements['ID_' + i] = struct.ProteinInfo.Elements.pop(i)
    for i in struct.ReactionInfo.Elements.keys():
        struct.ReactionInfo.Elements['ID_' + i] = struct.ReactionInfo.Elements.pop(i)
    for i in struct.MetaboliteConstraintsInfo.Elements.keys():
        struct.MetaboliteConstraintsInfo.Elements['ID_' +
                                                  i] = struct.MetaboliteConstraintsInfo.Elements.pop(i)
    for i in struct.EnzymeConstraintsInfo.Elements.keys():
        struct.EnzymeConstraintsInfo.Elements['ID_' +
                                              i] = struct.EnzymeConstraintsInfo.Elements.pop(i)
    for i in struct.ProcessConstraintsInfo.Elements.keys():
        struct.ProcessConstraintsInfo.Elements['ID_' +
                                               i] = struct.ProcessConstraintsInfo.Elements.pop(i)
    for i in struct.DensityConstraintsInfo.Elements.keys():
        struct.DensityConstraintsInfo.Elements['ID_' +
                                               i] = struct.DensityConstraintsInfo.Elements.pop(i)

    struct.ProcessInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.CompartmentInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.MetaboliteInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ModuleInfo.Elements.update({'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.EnzymeInfo.Elements.update({'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ProteinInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ReactionInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.DensityConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.ProcessConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.MetaboliteConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})
    struct.EnzymeConstraintsInfo.Elements.update(
        {'Description': 'abcdefghijklmnopqrstuvwxyz', 'Pictures': []})

    Block = {'ModelInformation': struct.GeneralInfo.JSONize(),
             'ModelStatistics': struct.ModelStatistics.JSONize(),
             'Process': struct.ProcessInfo.JSONize(),
             'Compartment': struct.CompartmentInfo.JSONize(),
             'Metabolite': struct.MetaboliteInfo.JSONize(),
             'Module': struct.ModuleInfo.JSONize(),
             'Enzyme': struct.EnzymeInfo.JSONize(),
             'Protein': struct.ProteinInfo.JSONize(),
             'Reaction': struct.ReactionInfo.JSONize(),
             'DensityConstraint': struct.DensityConstraintsInfo.JSONize(),
             'ProcessConstraint': struct.ProcessConstraintsInfo.JSONize(),
             'MetaboliteConstraint': struct.MetaboliteConstraintsInfo.JSONize(),
             'EnzymeConstraint': struct.EnzymeConstraintsInfo.JSONize()
             }

    return({'RBA_ModelData': {'StructuralInformation': Block}})


def generateOverview(StatsReactions, StatsMetabolites, StatsModules, StatsEnzymes, StatsProteins, StatsMacromolecules, StatsProcesses, StatsCompartments, StatsTargets, StatsConstraintsBiological, StatsConstraintsMathematical):
    out = {'Reactions Total': StatsReactions['ReactionsTotal'],
           'Reactions Unique': StatsReactions['ReactionsUnique'],
           'Reactions Spontaneous': StatsReactions['ReactionsSpontaneous'],
           'Reactions Enzymatic': StatsReactions['ReactionsEnzymatic'],
           'Reactions Internal': StatsReactions['ReactionsInternal'],
           'Reactions Exchange': StatsReactions['ReactionsExchange'],
           'Reactions CompartmentTransport': StatsReactions['ReactionsCompartmentTransport'],
           'Reactions Reversible': StatsReactions['ReactionsReversible'],
           'Reactions Irreversible': StatsReactions['ReactionsIrreversible'],
           'Metabolites Total': StatsMetabolites['MetabolitesTotal'],
           'Metabolites Internal': StatsMetabolites['MetabolitesInternal'],
           'Metabolites External': StatsMetabolites['MetabolitesExternal'],
           'Metabolites GrowthRelevant': StatsMetabolites['MetabolitesGrowthRelevant'],
           'Boundary Metabolites': StatsMetabolites['BoundaryMetabolites'],
           'Enzymes Total': StatsEnzymes['EnzymesTotal'],
           'Enzymes Unique': StatsEnzymes['EnzymesUnique'],
           'Proteins Total': StatsProteins['ProteinsTotal'],
           'RNAs Total': StatsMacromolecules['RNAsTotal'],
           'DNAs Total': StatsMacromolecules['DNAsTotal'],
           'Processes Total': StatsProcesses['ProcessesTotal'],
           'Modules Total': StatsModules['ModulesTotal'],
           'Compartments Total': StatsCompartments['CompartmentsTotal'],
           'Biological constraints metabolite': StatsConstraintsBiological['BioConstraintsMetabolite'],
           'Biological constraints capacity': StatsConstraintsBiological['BioConstraintsCapacity'],
           'Biological constraints process': StatsConstraintsBiological['BioConstraintsProcess'],
           'Biological constraints density': StatsConstraintsBiological['BioConstraintsDensity'],
           'Mathematical constraints variables': StatsConstraintsMathematical['MathConstraintsVariables'],
           'Mathematical constraints constraints': StatsConstraintsMathematical['MathConstraintsConstraints'],
           'Mathematical constraints equalities': StatsConstraintsMathematical['MathConstraintsEqualities'],
           'Mathematical constraints inequalities': StatsConstraintsMathematical['MathConstraintsInequalities']}
    out.update(StatsTargets)
    return(out)


def StatsConstraintsBiological(MetCs, CapCs, DenCs, EffCs):
    out = {'BioConstraintsMetabolite': len(MetCs.keys()),
           'BioConstraintsCapacity': len(CapCs.keys()),
           'BioConstraintsProcess': len(EffCs.keys()),
           'BioConstraintsDensity': len(DenCs.keys())}
    return(out)


def StatsConstraintsMathematical(MetCs, CapCs, DenCs, EffCs, Enzymes, Reactions, Processes):
    nVars = len(Enzymes.keys())+len(Reactions.keys())+len(Processes.keys())
    nConsts = len(MetCs.keys())+len(CapCs.keys())+len(DenCs.keys())+len(EffCs.keys())
    nEqC = 0
    nInC = 0
    for i in MetCs.keys():
        if MetCs[i]['Type'] == '<=':
            nInC += 1
        if MetCs[i]['Type'] == '=':
            nEqC += 1
    for i in CapCs.keys():
        if CapCs[i]['Type'] == '<=':
            nInC += 1
        if CapCs[i]['Type'] == '=':
            nEqC += 1
    for i in DenCs.keys():
        if DenCs[i]['Type'] == '<=':
            nInC += 1
        if DenCs[i]['Type'] == '=':
            nEqC += 1
    for i in EffCs.keys():
        if EffCs[i]['Type'] == '<=':
            nInC += 1
        if EffCs[i]['Type'] == '=':
            nEqC += 1
    out = {'MathConstraintsVariables': nVars,
           'MathConstraintsConstraints': nConsts,
           'MathConstraintsEqualities': nEqC,
           'MathConstraintsInequalities': nInC}
    return(out)


def findContainedProteins(comp, Prots):
    out = []
    for k in Prots.keys():
        if Prots[k]['Compartment'] == comp:
            out.append(k)
    return(out)


def findContainedMacromolecules(comp, Macromolecules):
    out = []
    for k in Macromolecules.keys():
        if Macromolecules[k]['Compartment'] == comp:
            out.append(k)
    return(out)


def findContainedEnzymes(ContainedProteins, Prots, Enzymes):
    rs = []
    enzs = []
    for k in ContainedProteins:
        rs = rs+Prots[k]['associatedReactions']
        enzs = enzs+Prots[k]['associatedEnzymes']
    out = {'containedEnzymes': list(numpy.unique(enzs)),
           'containedReactions': list(numpy.unique(rs))}
    return(out)


def findAssociatedEnzyme(Enzymes, protein):
    out1 = []
    out2 = []
    for e in Enzymes.keys():
        if protein in Enzymes[e]['Subunits'].keys():
            out1.append(e)
            out2.append(Enzymes[e]['Reaction'])
    out = {'Enz': out1,
           'Rx': out2}
    return(out)


def findIsozymes(ez, blocks, Reactions, rx):
    out = []
    twins = Reactions[rx]['Twins']
    if len(twins) > 0:
        for r in twins:
            if not type(r) == list:
                if not type(Reactions[r]['Enzyme']) == list:
                    out.append(Reactions[r]['Enzyme'])
    return(out)


def sortConstraints(matrix, model):
    metaboliteList = [model.metabolism.species._elements[m].id for m in range(
        len(model.metabolism.species._elements))]
    mets = {}
    capacity = {}
    efficiency = {}
    density = {}
    for j in range(len(matrix.row_names)):
        i = matrix.row_names[j]
        if i in metaboliteList:
            mets[i] = j
        if 'enzyme' in i:
            capacity[i] = j
        if i.startswith('P_'):
            efficiency[i] = j
        if '_density' in i:
            density[i] = j
    out = {'MetaboliteConsts': mets,
           'ProcessConsts': efficiency,
           'EnzymeConsts': capacity,
           'DensityConsts': density}
    return(out)


def JSON_Int64_compensation(o):
    if isinstance(o, numpy.int64):
        return int(o)
#    raise TypeError
