# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
import json
import urllib3
import pandas
import libsbml
from scipy.sparse import lil_matrix, hstack


# package imports
from rba.core.constraint_blocks import ConstraintBlocks
from rbatools.element_block import ElementBlock


class ReactionBlock(ElementBlock):
    """
    Class holding information on the reactions in the model.

   Attributes
   ----------
   Elements : dict
       Each model-enzyme is represented by a key.
       The values, holding information on each enzyme, are dicts with predefined keys:
           'ID' : reaction ID in model (type str)
           'Compartment_Machinery' : Localisation of enzyme-subunit, catalysing this reaction. (type list)
           'Name' : name according to BiGG (type str)
           'OtherIDs' : Other reaction names (eg. BiGG, KEGG) (type list)
           'Formula' : Reaction formula as string (type str)
           'Reactants' : Which metabolites does this reaction consume of and many (type dict)
           'Products' : Which metabolites does this reaction produce of and many (type dict)
           'Type' : Type of reaction ('normal' or 'transport') (type str)
           'Compartment_Species' : Location of the metabolites involved with this reaction (type list)
           'Enzyme' : Enzyme catalysing this reaction (type str)
           'Twins' : Isoreactions of this reactions (catalysed by iso-enzymes) (type list)
    """

    def fromFiles(self, model, Info, ReactionAnnotations, sbml, metaboliteBlock):
        """
        Derive reaction-info from RBA-model.

        Input
        -----
        RBA-model
        Dataframe, holding BiGG-reaction information.

        Returns
        -------
        Dictionary with reaction-info.

        """
        blocks = ConstraintBlocks(model)
        full_S = build_S(list(blocks.metabolism.external + blocks.metabolism.internal),
                         model.metabolism.reactions).toarray()
        self.Elements = {}
        self.BiGGids = []
        index = 0
        http = urllib3.PoolManager()
        if type(sbml) is not str:
            if type(sbml.model) is libsbml.Model:
                sbmlIDMap = [reaction.id for reaction in sbml.model.reactions]
        reconstruction = Info.loc['Reconstruction', 'Value']
        for i in range(len(blocks.metabolism.reactions)):
            Reactants = associatedReactants(i, blocks, model, full_S)
            Products = associatedProducts(i, blocks, model, full_S)
            Reversibility = checkReversibility(i, blocks)
            CompartmentInfo = findCompartment(
                i, blocks, Reactants, Products, Reversibility, metaboliteBlock)
            Twins = findTwinRxns(i, blocks)
            protoID = deriveProto_ID(i, blocks, Twins)
            IDdict = {'ProtoID': protoID}
            reactionName = ' '
            if type(sbml) is not str:
                if type(sbml.model) is libsbml.Model:
                    if protoID in sbmlIDMap:
                        IDdict.update(getReactionAnnotationsFromSBML(
                            sbmlIDMap.index(protoID), sbml))
                        reactionName = sbml.model.reactions[sbmlIDMap.index(protoID)].name
                    elif 'R_'+protoID in sbmlIDMap:
                        IDdict.update(getReactionAnnotationsFromSBML(
                            sbmlIDMap.index(str('R_'+protoID)), sbml))
                        reactionName = sbml.model.reactions[sbmlIDMap.index(str('R_'+protoID))].name
            if type(ReactionAnnotations) is pandas.core.frame.DataFrame:
                IDdict.update(readReactionAnnotations(i, ReactionAnnotations, blocks))
            self.BiGGids.append(protoID)
            index += 1
            self.Elements[blocks.metabolism.reactions[i]] = {'ID': blocks.metabolism.reactions[i],
                                                             'Compartment_Machinery': [],
                                                             'Name': reactionName,
                                                             'OtherIDs': IDdict,
                                                             'index': index,
                                                             'Formula': Reactants['rSide'] + ' <=> ' + Products['pSide'],
                                                             'Reactants': Reactants['reactants'],
                                                             'Products': Products['products'],
                                                             'Reversible': Reversibility['Reversible'],
                                                             'Type': CompartmentInfo['type'],
                                                             'Compartment_Species': CompartmentInfo['comp'],
                                                             'Enzyme': findAssociatedEnzyme(i, blocks),
                                                             'Twins': Twins}

    def overview(self):
        """
        Derive statistics on reactions.

        Returns
        -------
        Dictionary with general numbers on reactions.

        """
        nTot = len(list(self.Elements.keys()))
        nUnique = numpy.nan
        nSpont = 0
        nEnzyme = 0
        nInternalTransport = 0
        nExchange = 0
        nInternal = 0
        nRev = 0
        nIrrev = 0
        BiGGIDs = []
        for i in list(self.Elements.keys()):
            if len(self.Elements[i]['Enzyme']) > 0:
                nEnzyme += 1
            else:
                nSpont += 1
            if self.Elements[i]['Type'] == 'Normal':
                nInternal += 1
            if self.Elements[i]['Type'] == 'Transport (-->)':
                if len(self.Elements[i]['Reactants']) == 0 or len(self.Elements[i]['Products']) == 0:
                    nExchange += 1
                else:
                    nInternalTransport += 1

            if self.Elements[i]['Type'] == 'Transport (<==>)':
                if len(self.Elements[i]['Reactants']) == 0 or len(self.Elements[i]['Products']) == 0:
                    nExchange += 1
                else:
                    nInternalTransport += 1

            if self.Elements[i]['Reversible']:
                nRev += 1
            else:
                nIrrev += 1
#               BiGGIDs.append(self.Elements[i]['OtherIDs']['BiGG.Reaction'])
        nUnique = len(list(numpy.unique(self.BiGGids)))
        out = {'ReactionsTotal': nTot,
               'ReactionsUnique': nUnique,
               'ReactionsSpontaneous': nSpont,
               'ReactionsEnzymatic': nEnzyme,
               'ReactionsInternal': nInternal,
               'ReactionsExchange': nExchange,
               'ReactionsCompartmentTransport': nInternalTransport,
               'ReactionsReversible': nRev,
               'ReactionsIrreversible': nIrrev}
        return(out)


def deriveProto_ID(reaction, blocks, Twins):
    """
    Derive BiGG-ID from reactionID.
    Relies on the assumption that reactionID in RBA-model equals an 'R_', followed by the BiGG-ID.
    When the reaction has twins (duplicated reactions due to isozymes),
    the one without an appended '_2', '_3' ... needs to be found first.

    Returns
    -------
    String with derived BiGG-ID.
    """
    if len(Twins) > 0:  # Check if isozymic reactions exist
        # Find "original" reaction (shortest name)
        b = min(Twins+[blocks.metabolism.reactions[reaction]], key=len)
        out = b
        # Remove 'M_'-prefix
        if b.startswith('R_'):
            out = b[2:]
    else:
        out = blocks.metabolism.reactions[reaction]
        if out.startswith('R_'):
            out = out[2:]
    return(out)


def findReactionAnnotations(rx_name, http, reconstruction):
    #     if not reconstruction is 'GSMM':
    #        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' + reconstruction +'/reactions/' + rx_name)
    if not reconstruction is 'GSMM':
        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' +
                                reconstruction + '/reactions/' + rx_name)
    else:
        response = http.request('GET', 'http://bigg.ucsd.edu/api/v2/models/' +
                                'universal' + '/reactions/' + rx_name)
    try:
        x = json.loads(response.data.decode('utf-8'))
        out = {'BiGG': 'NA', 'KEGG': 'NA', 'BioCyc': 'NA', 'Name': ''}
        out['BiGG'] = rx_name
        out['Name'] = x['name']
        if 'KEGG Reaction' in list(x['database_links'].keys()):
            out['KEGG'] = x['database_links']['KEGG Reaction'][0]['id']
            out['BioCyc'] = x['database_links']['BioCyc'][0]['id']
        return(out)
    except:
        return({'BiGG': 'NA', 'KEGG': 'NA', 'BioCyc': 'NA', 'Name': ''})
#     else:
#        return('NA')


def getReactionAnnotationsFromSBML(index, sbml):
    out = {}
    for a in sbml.model.reactions[index].getAnnotationString().split('\n'):
        if 'rdf:resource="http://identifiers.org/' in a:
            annotationType = a.split(
                'rdf:resource="http://identifiers.org/')[1].split('"/>')[0].split('/')
            out.update({annotationType[0]: annotationType[1]})
    return(out)


def readReactionAnnotations(r, ReactionAnnotations, blocks):
    AnnotationKeys = list(ReactionAnnotations)
    AnnotationIDs = [numpy.nan]*len(AnnotationKeys)
    for i in AnnotationKeys:
        if blocks.metabolism.reactions[r] in list(ReactionAnnotations[i]):
            reaction = list(ReactionAnnotations[i]).index(blocks.metabolism.reactions[r])
            AnnotationIDs = list(ReactionAnnotations.loc[reaction])
    return(dict(zip(AnnotationKeys, AnnotationIDs)))


def findRxnName(rxnid, reactionsBiGG):
    """
    Retreive (descriptive) name of reaction from BiGG-file.

    Returns
    -------
    String with name (when found), otherwise empty.
    """
    if rxnid in list(reactionsBiGG.index):
        return(str(reactionsBiGG.loc[rxnid]['name']))
    else:
        return(str(''))


def associatedReactants(i, blocks, model, Sfull):
    """
    Derive information of reactant-side of reaction.

    Returns
    -------
    'reactants': Dictionary with reactants and stoichiometric factors.
    'rSide': String, representing reactant-side of reaction-formula.
    """
    rxn = model.metabolism.reactions.get_by_id(blocks.metabolism.reactions[i])
    if rxn is not None:
        reactants = {i.species: transform_to_int(i.stoichiometry)
                     for i in rxn.reactants._elements if i is not None}
        eq = ''
        for i in reactants.keys():
            if reactants[i] != 0:
                eq += str(reactants[i])+' '+str(i)+' '
    else:
        reactants = {}
        eq = ''
    return({'reactants': reactants, 'rSide': eq})


def associatedProducts(i, blocks, model, Sfull):
    """
    Derive information of product-side of reaction.

    Returns
    -------
    'products': Dictionary with products and stoichiometric factors.
    'pSide': String, representing product-side of reaction-formula.
    """
    rxn = model.metabolism.reactions.get_by_id(blocks.metabolism.reactions[i])
    if rxn is not None:
        products = {i.species: transform_to_int(i.stoichiometry)
                    for i in rxn.products._elements if i is not None}
        eq = ''
        for i in products.keys():
            if products[i] != 0:
                eq += str(products[i])+' '+str(i)+' '
    else:
        products = {}
        eq = ''
    return({'products': products, 'pSide': eq})


def checkReversibility(rx, blocks):
    """
    Information on default reaction flux-bounds and reversibility.

    Returns
    -------
    Dictionary with numerical values on flux-bounds and boolean for reversibility.
    """
    LB = blocks.metabolism._lb[rx].__dict__['value']
    UB = blocks.metabolism._ub[rx].__dict__['value']
    out = {'Reversible': True,
           'UB': UB,
           'LB': LB}
    if LB == 0:
        out['Reversible'] = False
    return(out)


def findCompartment(rx, blocks, aR, aP, rR, metaboliteBlock):
    """
    Derive information on compartment aspects of the reaction.

    Returns
    -------
    'type': String 'Transport','Exchange' or 'Normal' (kind of reaction).
    'comp': compartment of SBML-model. arrow between compartments when reaction is 'Transport'.
    """
    r = list(aR['reactants'].keys())
    if len(r) > 0:
        rComp = list(set([metaboliteBlock.Elements[rc]['Compartment'] for rc in r]))
    else:
        rComp = []

    p = list(aP['products'].keys())
    if len(p) > 0:
        pComp = list(set([metaboliteBlock.Elements[pc]['Compartment'] for pc in p]))
    else:
        pComp = []

    if set(rComp) == set(pComp):
        typ = 'Normal'
        comp = rComp
    elif len(list(set(list(rComp+pComp)))) == 1:
        comp = list(set(rComp+pComp))
        typ = 'Exchange'
    else:
        typ = 'Transport (-->)'
        comp = list(set(rComp+pComp))
        if rR['Reversible'] == 'True':
            typ = 'Transport (<==>)'
    out = {'type': typ,
           'comp': comp}
    return(out)


def findAssociatedEnzyme(rx, blocks):
    """
    Return enzyme species, associated with reaction.

    Returns
    -------
    String with enzymeID, empty string if reaction is spontaneous.
    """
    RN = blocks.metabolism.reactions[rx]
    if RN in blocks.enzymes.reaction_catalyzed:
        return(blocks.enzymes.ids[blocks.enzymes.reaction_catalyzed.index(RN)])
    else:
        return(str(''))


def findTwinRxns(rx, blocks):
    """
    Find Twin reactions (identical reactions, catalyzed by different (iso-)enzymes)

    Returns
    -------
    List of iso-reactions.
    """
    out = []
    if 'duplicate' in blocks.metabolism.reactions[rx]:
        for x in blocks.metabolism.reactions:
            if 'duplicate' in x:
                if blocks.metabolism.reactions[rx].rsplit('_duplicate')[0] == x.rsplit('_duplicate')[0]:
                    if blocks.metabolism.reactions[rx] is not x:
                        out.append(x)
            else:
                if blocks.metabolism.reactions[rx].rsplit('_duplicate')[0] == x:
                    if blocks.metabolism.reactions[rx] is not x:
                        out.append(x)
    else:
        for x in blocks.metabolism.reactions:
            if blocks.metabolism.reactions[rx]+'_duplicate' in x:
                out.append(x)
    return(out)


def build_S(metabolites, reactions):
    """
    Build stoichiometry matrix from metabolites and reactions.

    Parameters
    ----------
    metabolites:
        Metabolite identifiers (used to define row order).
    reactions: rba.xml.ListOfReactions
        Reaction data.

    Returns
    -------
    scipy.sparse.lil_matrix
        Stoichiometry matrix.

    """
    m_index = {m: i for i, m in enumerate(metabolites)}
    S = lil_matrix((len(metabolites), len(reactions)))
    for r_index, reaction in enumerate(reactions):
        for reactant in reaction.reactants:
            S[m_index[reactant.species], r_index] = -reactant.stoichiometry
        for product in reaction.products:
            S[m_index[product.species], r_index] = product.stoichiometry
    return S


def transform_to_int(number):
    if number % 1 == 0:
        return(int(number))
    else:
        return(number)
