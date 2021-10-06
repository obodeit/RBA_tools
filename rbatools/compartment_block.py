# python 2/3 compatibility
from __future__ import division, print_function

# package imports
from rbatools.element_block import ElementBlock


class CompartmentBlock(ElementBlock):
    """
    Class holding information on the compartments in the model.

    Attributes
    ----------
    Elements : dict
        Each model-compartment is represented by a key.
        The values, holding information on each compartment, are dicts with predefined keys:
            'ID' : compartment ID in model (type str)
            'associatedProteins' : proteins localised to compartment (type list)
            'associatedEnzymes' : enzyme located in compartment (type list)
            'associatedReactions' : metabolic reactions located in compartment (type list)
    """

    def fromFiles(self, model, Info):

        Compartments = getCompartmentList(model)
        self.Elements = {}
        index = 0
        for i in Compartments:
            index += 1
            self.Elements[i] = {'ID': i,
                                'associatedProteins': [],
                                'index': index,
                                'associatedReactions': [],
                                'associatedEnzymes': []}

    def overview(self):
        """
        Derive statistics on compartments.

        Returns
        -------
        Dictionary with general numbers on compartments.

        """
        nT = len(self.Elements.keys())
        out = {'CompartmentsTotal': nT}
        return(out)


def getCompartmentList(model):
    out = []
    for c in range(len(model.metabolism.compartments._elements)):
        out.append(model.metabolism.compartments._elements[c].id)
    return(out)
