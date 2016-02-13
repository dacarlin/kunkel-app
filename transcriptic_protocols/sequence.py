import collections
from autoprotocol import UserError
from modules.utils import *


def sequence(protocol, params):

    seq_primers = collections.OrderedDict({})
    for well_set in params['seq_set']:
        for primer in well_set['seq_primers']:
            if primer in seq_primers.keys():
                seq_primers[primer].extend(well_set['growth_wells'].wells)
            else:
                seq_primers[primer] = []
                seq_primers[primer].extend(well_set['growth_wells'].wells)

    seq_plates = []
    for i, primer in enumerate(seq_primers.keys()):
        if len(seq_primers[primer]) >= 97:
            raise UserError('You can only sequence up to 96 wells per primer,'
                            ' please submit separate sequencing protocols.')
        seq_plates.append(protocol.ref('seq_plate_%s_%s' % (primer.container.name, printdatetime(time=False)), cont_type='96-pcr', storage='cold_4'))
        seq_wells = seq_plates[i].all_wells(columnwise=True)
        for j, well in enumerate(seq_primers[primer]):
            protocol.transfer(well, seq_wells[j], "30:microliter")
            if well.name:
                seq_wells[j].set_name(well.name)
            else:
                seq_wells[j].set_name('well_%s_%s' % (well.humanize(), str(primer.container.name)))

    for k, primer in enumerate(seq_primers.keys()):
        dataref = "seq_primer:_%s" % (primer.container.name,)
        protocol.sangerseq(seq_plates[k], seq_plates[k].wells_from(0, len(seq_primers[primer]), columnwise=True).indices(), dataref, type="rca",
                               primer=primer.container)


if __name__ == '__main__':
    from autoprotocol.harness import run
    run(sequence, 'Sequence')
