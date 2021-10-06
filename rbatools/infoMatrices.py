# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
import pandas
import json

class InfoMatrices(object):
     """
     Class holding information on the compartments in the model.

     Attributes
     ----------
     Reaction_Reaction : pandas.DataFrame
          Mapping of model-reactions to unique BiGG IDs
          (compensation for reaction-duplication for Isozymes).
          rows: unique BiGG-reactions ; cols: model-reactions
     Protein_Enzyme : pandas.DataFrame
          Conversion of Enzyme levels to levels of proteins, constituting them.
          (Gives Proteome of Metabolic Proteins, needs to be added to Process_Protein to obtain full proteome)
          rows: proteins ; cols: model-enzymes
     ProteinWeight : pandas.DataFrame
          Molecular weight and number of amino-acids for each protein
          rows: proteins ; cols: AAnum and weight
     Protein_ProcessMachinery : pandas.DataFrame
          Conversion of Process levels to levels of proteins, constituting their machineries.
          (Gives Proteome of Process Proteins, needs to be added to Enzyme_Protein to obtain full proteome)
          rows: proteins ; cols: model-processes
     Process_Protein : pandas.DataFrame
          Indication on how much each protein requires of each process
          rows: processes ; cols: protein
     Compartment_Protein : pandas.DataFrame
          Logical matrix on which protein is located in which compartment
          rows: compartments ; cols: proteins
     S : pandas.DataFrame
          Model stoichiometry-matrix
          rows: metabolites ; cols: reactions

     """
     def __init__(self,struct):
          BiGGids=[]
          for rx in list(struct.ReactionInfo.Elements.keys()):
                Bid=struct.ReactionInfo.Elements[rx]['OtherIDs']['ProtoID']
                if Bid is not '':
                     BiGGids.append(Bid)

          self.Reaction_Reaction = makeReaction_Reaction(struct,BiGGids)
          self.Protein_Enzyme = makeProtein_Enzyme(struct)['Matrix']
          self.ProteinWeight = makeProtein_Enzyme(struct)['Weight']
          self.Protein_ProcessMachinery = makeProcessMachinery_Protein(struct)
          self.Process_Protein = makeProcessRequirements_Protein(struct)
          self.Compartment_Protein = makeCompartment_Protein(struct)
          self.S = make_S(struct)

def make_S(input):
     Mets=input.MetaboliteInfo.toDataFrame()
     InternalMets=sorted(list(Mets[Mets['Type'].apply(json.loads)=='internal']['ID'].apply(json.loads)))
     ExternalMets=sorted(list(Mets[Mets['Type'].apply(json.loads)=='external']['ID'].apply(json.loads)))
     PrecursorMets=sorted(list(Mets[Mets['Type'].apply(json.loads)=='precursor']['ID'].apply(json.loads)))
     Metabolites=ExternalMets + InternalMets + PrecursorMets
     Rxns=input.ReactionInfo.toDataFrame()
     Reactions=sorted(Rxns.index.tolist())
     S=pandas.DataFrame(numpy.zeros((len(Metabolites),len(Reactions))),index=Metabolites,columns=Reactions)
     for rx in Reactions:
          for j in list(json.loads(Rxns.loc[rx]['Reactants']).keys()):
               S.loc[j,rx]=-json.loads(Rxns.loc[rx]['Reactants'])[j]
          for j in list(json.loads(Rxns.loc[rx]['Products']).keys()):
               S.loc[j,rx]=json.loads(Rxns.loc[rx]['Products'])[j]
     return(S)

def makeReaction_Reaction(struct,BiGGids):
     R_Matrix=pandas.DataFrame(numpy.zeros((len(numpy.unique(BiGGids)),len(list(struct.ReactionInfo.Elements.keys())))),columns=list(struct.ReactionInfo.Elements.keys()),index=numpy.unique(BiGGids))
     for rx in list(struct.ReactionInfo.Elements.keys()):
          R_Matrix.loc[struct.ReactionInfo.Elements[rx]['OtherIDs']['ProtoID'],rx]=1
     return(R_Matrix)

def makeProtein_Enzyme(struct):
     P_Matrix=pandas.DataFrame(numpy.zeros((len(numpy.unique(list(struct.ProteinInfo.Elements.keys()))),len(numpy.unique(list(struct.EnzymeInfo.Elements.keys()))))),columns=numpy.unique(list(struct.EnzymeInfo.Elements.keys())),index=numpy.unique(list(struct.ProteinInfo.Elements.keys())))
     Pweight=pandas.DataFrame(numpy.zeros((len(numpy.unique(list(struct.ProteinInfo.Elements.keys()))),2)),columns=['AAlength','MolecMass'],index=numpy.unique(list(struct.ProteinInfo.Elements.keys())))
     for i in numpy.unique(list(struct.ProteinInfo.Elements.keys())):
          Pweight.loc[i,'AAlength']=struct.ProteinInfo.Elements[i]['AAnumber']
          Pweight.loc[i,'MolecMass']=struct.ProteinInfo.Elements[i]['Weight']
          if len(struct.ProteinInfo.Elements[i]['associatedEnzymes']) >0:
               for j in struct.ProteinInfo.Elements[i]['associatedEnzymes']:
                    P_Matrix.loc[i,j]=struct.EnzymeInfo.Elements[j]['Subunits'][i]['StochFac']
     return({'Matrix': P_Matrix , 'Weight': Pweight})

def makeProcessRequirements_Protein(struct):
     PM_Matrix=pandas.DataFrame(numpy.zeros((2,len(numpy.unique(list(struct.ProteinInfo.Elements.keys()))))),columns=numpy.unique(list(struct.ProteinInfo.Elements.keys())),index=['P_TA','P_CHP'])
     for i in numpy.unique(list(struct.ProteinInfo.Elements.keys())):
          t=0
          f=0
          if 'Translation' in list(struct.ProteinInfo.Elements[i]['ProcessRequirements'].keys()):
              t=struct.ProteinInfo.Elements[i]['ProcessRequirements']['Translation']
          if 'Folding' in list(struct.ProteinInfo.Elements[i]['ProcessRequirements'].keys()):
              f=struct.ProteinInfo.Elements[i]['ProcessRequirements']['Folding']
          PM_Matrix.loc['P_TA',i]=t
          PM_Matrix.loc['P_CHP',i]=f
     return(PM_Matrix)

def makeProcessMachinery_Protein(struct):
     M_Matrix=pandas.DataFrame(numpy.zeros((len(numpy.unique(list(struct.ProteinInfo.Elements.keys()))),2)),columns=['P_TA','P_CHP'],index=numpy.unique(list(struct.ProteinInfo.Elements.keys())))
     for i in numpy.unique(list(struct.ProteinInfo.Elements.keys())):
          if i in list(struct.ProcessInfo.Elements['Translation']['Composition'].keys()):
               M_Matrix.loc[i,'P_TA']=struct.ProcessInfo.Elements['Translation']['Composition'][i]
          if i in list(struct.ProcessInfo.Elements['Folding']['Composition'].keys()):
               M_Matrix.loc[i,'P_CHP']=struct.ProcessInfo.Elements['Folding']['Composition'][i]
     return(M_Matrix)

def makeCompartment_Protein(struct):
     CP_Matrix=pandas.DataFrame(numpy.zeros((len(numpy.unique(list(struct.CompartmentInfo.Elements.keys()))),len(numpy.unique(list(struct.ProteinInfo.Elements.keys()))))),columns=numpy.unique(list(struct.ProteinInfo.Elements.keys())),index=numpy.unique(list(struct.CompartmentInfo.Elements.keys())))
     for i in numpy.unique(list(struct.CompartmentInfo.Elements.keys())):
          CP_Matrix.loc[i,struct.CompartmentInfo.Elements[i]['associatedProteins']]=1
     return(CP_Matrix)
