# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
import libsbml
# package imports
from rba.core.constraint_blocks import ConstraintBlocks
from rbatools.element_block import ElementBlock


class ModuleBlock(ElementBlock):
    """
    Class holding information on the modules in the model.

  Attributes
  ----------
  Elements : dict
      Each model-enzyme is represented by a key.
      The values, holding information on each enzyme, are dicts with predefined keys:
          'ID' : So far dummy (type str)
          'Name' : So far dummy (type str)
    """

    def fromFiles(self, model, sbml):
        self.Elements = {}
        if type(sbml) is not str:
            if type(sbml.model) is libsbml.Model:
                GroupsOfSBML = sbml.model.getPlugin('groups')
                if GroupsOfSBML is not None:
                    if len(GroupsOfSBML.getListOfGroups()) > 0:
                        for group in GroupsOfSBML.getListOfGroups():
                            groupToAdd = {}
                            groupToAdd['ID'] = group.getIdAttribute()
                            groupToAdd['Name'] = group.getName()
                            groupToAdd['Members'] = [member.getIdRef()
                                                     for member in group.getListOfMembers()]
                            self.Elements.update({groupToAdd['ID']: groupToAdd})

    def overview(self):
        """
        Derive statistics on modules.

        Returns
        -------
        Dictionary with general numbers on modules.

        """
        nT = len(self.Elements.keys())
        out = {'ModulesTotal': nT}
        return(out)
