# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import numpy
# package imports
from rbatools.element_block import ElementBlock


class TargetBlock(ElementBlock):
    """
    Class holding information on the targets in the model.

   Attributes
   ----------
   Elements : dict
       Each model-target is represented by a key.
       The values, holding information on each target, are dicts with predefined keys:
           'ID' : Target ID in model (type str)
           'Group' : Class of targets (type str):
                Eg:
                'macrocomponent_production','maintenance_atp_target',
                'metabolite_production','replication_targets',
                'rna_degradation' or 'translation_targets'
           'Type' : Type of targets (type str):
                Eg:
                'metabolite_production', 'degradation_fluxes',
                'concentrations'  or 'reaction_fluxes'
           'TargetSpecies' : For which entity the target is defined (type str)
           'TargetValue' : Model-parameter ID which defines target value (type str)
        """

    def fromFiles(self, model):
        self.Elements = {}
        self.GroupList = []
        nGroups = len(list(model.targets.target_groups._elements))
        for i in range(nGroups):
            typeList = list(model.targets.target_groups._elements[i].__dict__.keys())
            typeList.remove('id')
            for j in typeList:
                for k in model.targets.target_groups._elements[i].__dict__[j].__dict__['_elements']:
                    self.GroupList.append(
                        str(model.targets.target_groups._elements[i].__dict__['id']))
                    spec = 'Other'
                    if 'species' in list(k.__dict__.keys()):
                        spec = str(k.__dict__['species'])
                    if 'reaction' in list(k.__dict__.keys()):
                        spec = str(k.__dict__['reaction'])
                    if str(j) == 'degradation_fluxes':
                        targettype_specifier = 'degradation'
                    elif str(j) == 'production_fluxes':
                        targettype_specifier = 'production'
                    elif str(j) == 'reaction_fluxes':
                        targettype_specifier = 'flux'
                    elif str(j) == 'concentrations':
                        targettype_specifier = 'concentration'
                    else:
                        targettype_specifier = ''
                    Edict = {'ID': 'Target_'+targettype_specifier+'_'+spec,
                             'Group': str(model.targets.target_groups._elements[i].__dict__['id']),
                             'Type': str(j),
                             'TargetEntity': spec,
                             'TargetValue': str(k.__dict__['value'])}
                    self.Elements.update({Edict['ID']: Edict})

    def overview(self):
        """
        Derive statistics on targets.

        Returns
        -------
        Dictionary with general numbers on targets.

        """
        out = {}
        for i in numpy.unique(self.GroupList):
            nI = len(list(numpy.where(numpy.array(self.GroupList) == i)[0]))
            out.update({'Targets_'+i: nI})
        return(out)


# model.targets.target_groups._elements --> list of different target groups ("translation_targets", "transcription_targets" ...)


# model.targets.target_groups._elements[0].__dict__['id'] --> Group
# list(model.targets.target_groups._elements[0].__dict__.keys()).remove('id') --> types
