from autoprotocol.protocol import Protocol
from autoprotocol.instruction import Instruction
from assemble import assemble
from sequence import sequence
from transform import transform
from plasmidprep import plasmidprep
from autoprotocol.container import WellGroup
from mutant import Mutant
from collections import OrderedDict
from autoprotocol import UserError
from modules.utils import printdatetime

def scale_default(length, scale, label):
    ok = True
    if scale == '25nm':
        ok = True if (length >= 15 and length <= 60) else False
    elif scale == '100nm':
        ok = True if (length >= 10 and length <= 90) else False
    elif scale == '250nm':
        ok = True if (length >= 5 and length <= 100) else False
    elif scale == '1um':
        ok = True if (length >= 5 and length <= 100) else False
    else:
        ok = False
    if not ok:
        raise UserError("The specified oligo '%s' is %s base pairs long."
                           " This sequence length is invalid for the scale "
                           "of synthesis chosen (%s)." % (label, length, scale))


def kunkel_full(protocol, params):
    growth_media = params["construct_setup"]['growth_media']
    num_colonies = params["construct_setup"]['num_colonies']
    ssDNA = params["construct_setup"]['ssDNA']
    mutant_constructs = []

    # make mutant objects for accessibility
    construct_collect = {}
    for csv_row in params["construct_setup"]['mutant_upload']:
        if csv_row["mutant_label"] not in construct_collect.keys():
            construct_collect[csv_row["mutant_label"]] = []
            construct_collect[csv_row["mutant_label"]].append(
                {
                "sequence": csv_row["sequence"],
                "purification": csv_row["purification"],
                "scale": csv_row["scale"],
                "oligo_label": csv_row["oligo_label"]

            })
        else:
            construct_collect[csv_row["mutant_label"]].append(
                {
                "sequence": csv_row["sequence"],
                "purification": csv_row["purification"],
                "scale": csv_row["scale"],
                "oligo_label": csv_row["oligo_label"]

            }
                )

    oligo_collect = {}
    for row in params["construct_setup"]["mutant_upload"]:
        if (row["sequence"] not in oligo_collect.keys() and row["oligo_label"] in protocol.refs.keys()):
            raise RuntimeError("You cannot specify two different "
                   "oligos to be synthesized with the "
                   "same name %s" % row['oligo_label'])
        elif row["sequence"] not in oligo_collect.keys():
            oligo_collect[row["sequence"]] = {
                "sequence": row["sequence"],
                "purification": row["purification"],
                "scale": row["scale"],
                "destination": protocol.ref(row["oligo_label"], None, "micro-2.0", storage="cold_4").well(0)
            }


    for mut in construct_collect.keys():
        mut_oligos = [o for o in construct_collect[mut]]
        mutant = Mutant(mut)
        for oligo in mut_oligos:
            mutant.add_oligos(oligo_collect[oligo["sequence"]]["destination"])
        mutant_constructs.append(mutant)


    oligos_to_synthesize = []
    for o in oligo_collect.keys():
        scale_default(len(oligo_collect[o]["sequence"]), oligo_collect[o]["scale"], oligo_collect[o]["destination"].container.name)
        oligos_to_synthesize.append(oligo_collect[o])
    protocol.oligosynthesize(oligos_to_synthesize)

    assemble_params = {
        'ssDNA': ssDNA,
        'constructs': [{
            'mutant_name': mu.name,
            'oligos': mu.oligos} for mu in mutant_constructs],
        'mutant_objs': mutant_constructs

    }

    annealing_plate = assemble(protocol, assemble_params)
    protocol.unseal(annealing_plate)

    transform_params = {
        'num_colonies': num_colonies,
        'growth_media': growth_media,
        'constructs': [mu.anneal_well for mu in mutant_constructs],
        'mutant_objs': mutant_constructs

    }

    growth_plate = transform(protocol, transform_params)

    seq_primers = []

    for seq_primer in params["sequencing"]:
        if seq_primer["seq_choice"] != "No sequencing.":
        # make temp container with name of stock primer
            primer = protocol.ref(seq_primer["seq_choice"], None, "micro-1.5",
                                  discard=True).well(0)
            primer.set_name(seq_primer["seq_choice"])
            primer_vol = "1:microliter"
            seq_primers.append(primer)

    sequence_params = {
        'seq_set': [{
            'growth_wells': WellGroup([w for w in mu.growth_wells]),
            'seq_primers': seq_primers} for mu in mutant_constructs]

    }
    if seq_primers:
        protocol.uncover(growth_plate)
        sequence(protocol, sequence_params)
        protocol.cover(growth_plate, lid="low_evaporation")

    if params["other_processing"]["other_processing"] != "No processing.":
        protocol.uncover(growth_plate)
        if params["other_processing"]["other_processing"] == "Miniprep":
            mini_samples = []
            for mu in mutant_constructs:
                for w in mu.growth_wells:
                    mini_samples.append({"sample": w, "name": w.name})

            miniprep_params = {
                        "type": "Miniprep",
                        "media": growth_media,
                        'samples': mini_samples,
                        "growth_plate": growth_plate
                        }

            plasmidprep(protocol, miniprep_params)

        if params["other_processing"]["other_processing"] == "Return Colonies":
            return_plate = protocol.ref("return_plate_%s" % printdatetime(time=False), cont_type='96-pcr', storage='cold_4')
            for mut in  mutant_constructs:
                for g_well in mut.growth_wells:
                    protocol.transfer(g_well, return_plate.well(g_well.index), "30:microliter")
                    return_plate.well(g_well.index).set_name(g_well.name)

            protocol.seal(return_plate)
            protocol.cover(growth_plate, lid="low_evaporation")

        if params["other_processing"]["other_processing"] == "Return Colonies Glycerol":
            return_plate = protocol.ref("return_plate_glycerol_%s" % printdatetime(time=False), cont_type='96-pcr', storage='cold_4')
            for mut in  mutant_constructs:
                for g_well in mut.growth_wells:
                    protocol.transfer(g_well, return_plate.well(g_well.index), "30:microliter")
                    return_plate.well(g_well.index).set_name(g_well.name)
            for mut in  mutant_constructs:
                for g_well in mut.growth_wells:
                    protocol.provision("rs17rrhqpsxyh2", return_plate.well(g_well.index), "30:microliter")

            protocol.seal(return_plate)
            protocol.cover(growth_plate, lid="low_evaporation")


if __name__ == '__main__':
    from autoprotocol.harness import run
    run(kunkel_full, 'KunkelSiegel')
