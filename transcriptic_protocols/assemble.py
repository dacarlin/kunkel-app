import math
from modules.utils import provision_to_tube, printdatetime, det_new_group, thermocycle_ramp
from autoprotocol.pipette_tools import aspirate_source, depth


def assemble(protocol, params):

    def provision_reagents(reagents, dest, num_rxts, mm_mult=1.3, num_rxts_plus=3.0):
        for reagent in reagents.values():
            protocol.provision(reagent['resource_id'], dest, "%s:microliter" % ((num_rxts + num_rxts_plus) * reagent['reagent_ratio'] * mm_mult))

    def transfer_kwargs(pre_buffer, one_tip=False, one_source=False):
        kwargs = {"one_tip": one_tip,
              "one_source": one_source,
              "pre_buffer": "%s:microliter" % pre_buffer,
              "blowout_buffer": True}
        return(kwargs)

    # general parameters
    ssDNA = params['ssDNA']
    constructs = [construct['oligos'] for construct in params['constructs']]
    num_constructs = len(constructs)
    flattened = [val for oligo in constructs for val in oligo]
    oligos = list({v: v for v in flattened}.values())
    num_oligos = len(oligos)
    mm_mult = 1.3
    num_rxts_plus = 3

    # provision water
    water = [provision_to_tube(protocol, "water%s" % (i + 1), "micro-2.0", "rs17gmh5wafm5p", 1900)
             for i in range(int(math.ceil(num_constructs/float(9.0))))
            ]

    for w in water:
        w.set_volume('1800:microliter')

    # provision atp for entire protocol
    atp_reagents = {'atp': {"resource_id": 'rs16pccshb6cb4',
                    'reagent_ratio': 0.1},
                    'water': {"resource_id": 'rs17gmh5wafm5p', 'reagent_ratio': 0.9}}

    atp_rxts = (num_oligos + num_constructs)
    atp = protocol.ref("atp_10mM", cont_type='micro-1.5', discard=True).well(0)
    provision_reagents(atp_reagents, atp, atp_rxts, mm_mult, num_rxts_plus=6)

    # provision kinase mix
    kinase_mix = []
    for i in range(int(math.ceil(num_oligos/60.0))):
        kinase_mix.append(protocol.ref("kinase_mix-%s" % (i + 1), None, "micro-1.5", discard=True).well(0))

    reagents = {'pnkbuffer': {"resource_id": 'rs16pc9rd5sg5d', "reagent_ratio": 3},
                'water': {"resource_id": 'rs17gmh5wafm5p', "reagent_ratio": 18},
                'pnk': {"resource_id": 'rs16pc9rd5hsf6', "reagent_ratio": 1}}

    provision_reagents(reagents, kinase_mix, num_oligos, mm_mult, num_rxts_plus)
    protocol.transfer(atp, kinase_mix, "%s:microliter" % ((num_oligos + num_rxts_plus) * 1 * mm_mult), new_group=True)

    #   kinase oligos
    kinase_oligo_plate = protocol.ref("kinase_oligo_plate_%s" % printdatetime(time=False), None, "96-pcr",
                                      discard=True)
    wells_to_kinase = kinase_oligo_plate.wells_from("A1", num_oligos)
    protocol.transfer(kinase_mix,
                      wells_to_kinase,
                      "23:microliter",
                      **transfer_kwargs(15, True, True))

    for i, oligo in enumerate(oligos):
        protocol.transfer(oligo,
                              wells_to_kinase[i],
                              "7:microliter",
                              mix_after=False,
                              new_group=det_new_group(i),
                              aspirate_source=aspirate_source(depth=depth("ll_following",
                                                              lld="pressure",
                                                              distance="0.0:meter")),
                              **transfer_kwargs(10))

    protocol.seal(kinase_oligo_plate)

    protocol.thermocycle(kinase_oligo_plate,
                         [{"cycles": 1, "steps": [
                             {"temperature": "37:celsius",
                              "duration": "60:minute"},
                             ]}
                          ], volume="30:microliter")

    # dilute oligos
    protocol.unseal(kinase_oligo_plate)

    diluted_oligo_plate = protocol.ref("dilute_oligo_plate", None, "96-flat", discard=True)
    diluted_oligo_wells = diluted_oligo_plate.wells_from(0, num_constructs)

    protocol.transfer(water,
                      diluted_oligo_wells,
                      "200:microliter",
                      disposal_vol="0:microliter",
                      **transfer_kwargs(40, True, True))

    for i, m in enumerate(constructs):
        for kin_oligo in m:
            index = list(oligos).index(kin_oligo)
            protocol.transfer(kinase_oligo_plate.well(index), diluted_oligo_plate.well(i),
                              "2:microliter",
                              mix_after=False,
                              mix_vol="5:microliter")
        diluted_oligo_plate.well(i).set_name(params['constructs'][i]['mutant_name'])

    protocol.cover(diluted_oligo_plate)
    protocol.spin(diluted_oligo_plate, "250:g", "2:minute")

    # make ssDNA_mastermix
    mm_mult_ssDNA = 1.5

    mix_plate = protocol.ref("mix_plate", None, "96-pcr", discard=True)
    ssDNA_mix = mix_plate.well(0)
    protocol.transfer(ssDNA,
                      ssDNA_mix,
                      "%s:microliter" % ((num_constructs + num_rxts_plus) * 2.0 * mm_mult_ssDNA),
                      **transfer_kwargs((num_constructs + 1) * 1))
    protocol.provision('rs17sh5rzz79ct', ssDNA_mix, "%s:microliter" % ((num_constructs + num_rxts_plus) * 0.2 * mm_mult_ssDNA))

    # anneal oligos
    protocol.uncover(diluted_oligo_plate)
    annealing_plate = protocol.ref("annealing_oligo_plate", None, "384-pcr", storage="cold_20")
    anneal_wells = annealing_plate.wells_from(0, num_constructs)
    protocol.transfer(ssDNA_mix,
                      anneal_wells.wells,
                      "2.2:microliter",
                      dispense_speed="50:microliter/second",
                      **transfer_kwargs(7, True, True))

    for dil_oligo, reaction in zip(diluted_oligo_wells.wells, anneal_wells.wells):
        protocol.transfer(dil_oligo,
                          reaction,
                          "2:microliter",
                          aspirate_source=aspirate_source(depth("ll_bottom", distance=".001:meter")),
                          mix_after=True,
                          mix_vol="2:microliter",
                          flowrate="50:microliter/second",
                          repetitions=2,
                          new_group=det_new_group(i),
                          **transfer_kwargs(5))
        reaction.set_name(dil_oligo.name)

    protocol.seal(annealing_plate)
    protocol.spin(annealing_plate, "250:g", "2:minute")
    protocol.thermocycle(annealing_plate, [{
        "cycles": 1,
        "steps": thermocycle_ramp("95:celsius", "25:celsius", "60:minute", "4:minute")
        }],
        volume="5:microliter",
        dataref=None,
        dyes=None)

    # polymerize
    protocol.unseal(annealing_plate)
    polymerize_MM = mix_plate.well(12)
    reagents = {"buffer": {"resource_id": 'rs17sh5rzz79ct', "reagent_ratio": 0.6},
                "t4ligase": {"resource_id": 'rs16pc8krr6ag7', "reagent_ratio": 0.4},
                "t7polymerase": {"resource_id": 'rs16pca2urcz74', "reagent_ratio": 0.4},
                "dntp": {"resource_id": 'rs16pcb542c5rd', "reagent_ratio": 0.4}
                }

    provision_reagents(reagents, polymerize_MM, num_constructs, mm_mult, num_rxts_plus)
    protocol.transfer(atp, polymerize_MM, "%s:microliter" % ((num_constructs + num_rxts_plus) * 0.4 * mm_mult), new_group=True)

    for reaction in anneal_wells.wells:
        protocol.transfer(polymerize_MM,
                          reaction,
                          "2.2:microliter",
                          mix_after=False,
                          **transfer_kwargs(10))
        if 'mutant_objs' in params.keys():
            mut = next(m for m in params['mutant_objs'] if m.name == reaction.name)
            mut.anneal_well = reaction

    protocol.seal(annealing_plate)
    protocol.incubate(annealing_plate, "ambient", "1.5:hour")

    # pass plate back for unsealing
    return annealing_plate


if __name__ == '__main__':
    from autoprotocol.harness import run
    run(assemble, "Assemble")
