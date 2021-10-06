from __future__ import division, print_function
import pandas
import json
from datetime import datetime


class RBA_LogBook(object):
    def __init__(self, type):
        self.type = type
        self.History = pandas.DataFrame(columns=['type', 'changes'])

    def addEntry(self, entry):
        intermediateEntry = pandas.DataFrame(columns=['type', 'changes'])
        tt = datetime.now().strftime("%Y/%m/%d, %H:%M:%S.%f")
        intermediateEntry.loc[tt, 'type'] = self.type
        intermediateEntry.loc[tt, 'changes'] = entry
        self.History = pandas.concat([self.History, intermediateEntry], axis=0)

    def mergeHistories(self, logbook):
        self.History = pandas.concat([self.History, logbook.History], axis=0).sort_index()
