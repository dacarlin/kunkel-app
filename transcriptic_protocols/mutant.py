from autoprotocol import UserError

class Mutant(object):

        unique_mutants = {}

        def __init__(self, name, oligos=None, seq_primers=None):
            if name in Mutant.unique_mutants:
                raise UserError("You must specify unique mutant names for"
                                " your constructs. %s has already been used." % name)
            self.oligos = oligos or []
            self.name = name
            self.seq_primers = seq_primers or []
            self.growth_wells = []
            Mutant.unique_mutants[self.name] = self

        def add_seq_primers(self, seq_primers):
            if type(seq_primers) is list:
                self.seq_primers.extend(seq_primers)
            else:
                self.seq_primers.append(seq_primers)

        def add_oligos(self, oligos):
            if type(oligos) is list:
                self.oligos.extend(oligos)
            else:
                self.oligos.append(oligos)

        def add_growth_wells(self, growth_wells):
            if type(growth_wells) is list:
                self.growth_wells.extend(growth_wells)
            else:
                self.growth_wells.append(growth_wells)
