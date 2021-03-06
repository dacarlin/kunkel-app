import math
from autoprotocol import UserError
from modules.utils import *


def transform(protocol, params):

    # general parameters
    constructs = params['constructs']
    num_constructs = len(constructs)
    plates = list(set([construct.container for construct in constructs]))
    if len(plates) != 1:
        raise UserError('You can only transform aliquots from one common container.')

    # **** need to be able to check if plate is sealed to add run-chaining ****
    num_colonies = params['num_colonies']
    mm_mult = 1.3

    # provision SOC medium
    soc_vol = num_constructs * 50 * mm_mult
    num_soc = int(math.ceil(float(soc_vol)/1800))
    soc_medium = [provision_to_tube(protocol, "soc_medium_%d" % (i+1),
                                    "micro-2.0", "rs17tpdy56hfar", soc_vol/num_soc).set_volume(Unit(0.9 * (soc_vol/num_soc), "microliter"))
                  for i in range(0, num_soc)]

    transformation_plate = protocol.ref("transformation_plate", None, "96-pcr", discard=True)
    protocol.incubate(transformation_plate, "cold_20", "10:minute")

    transformation_wells = transformation_plate.wells_from(0, num_constructs)
    for i in range(num_constructs):
        protocol.provision("rs16pbjc4r7vvz", transformation_wells[i], "50:microliter")

    for i, well in enumerate(constructs):
        protocol.transfer(well, transformation_wells[i], "2.0:microliter",
                          dispense_speed="10:microliter/second",
                          mix_after=False,
                          new_group=det_new_group(i))
        if well.name:
            transformation_wells[i].set_name(well.name)
        else:
            transformation_wells[i].set_name('construct_%s' % (i+1))

    # NEED to confirm second de-seal is working OR move to cover/uncover 96-flat
    protocol.seal(transformation_plate)
    protocol.incubate(transformation_plate, "cold_4", "20:minute", shaking=False, co2=0)
    protocol.unseal(transformation_plate)
    protocol.transfer(soc_medium, transformation_wells, "50:microliter", one_source=True)
    protocol.seal(transformation_plate)
    protocol.incubate(transformation_plate, "warm_37", "10:minute", shaking=True)
    protocol.unseal(transformation_plate)

    agar_plates = []
    agar_wells = WellGroup([])
    for well in range(0, len(transformation_wells), 6):
        agar_plate = ref_kit_container(protocol,
                                       "agar-%s_%s" % (len(agar_plates), printdatetime(time=False)),
                                       "6-flat",
                                       return_agar_plates(6)[params['growth_media']],
                                       discard=False, store='cold_4')
        agar_plates.append(agar_plate)
        for i, w in enumerate(transformation_wells[well:well + 6]):
            protocol.spread(w, agar_plate.well(i), "100:microliter")
            agar_wells.append(agar_plate.well(i).set_name(w.name))

    for agar_p in agar_plates:
        protocol.incubate(agar_p, "warm_37", "16:hour")

    growth_plate = protocol.ref("growth_plate_kunkel_%s" % printdatetime(time=False), None, "96-flat", storage='cold_4')
    total_cols = num_colonies * num_constructs
    if total_cols > 96:
        raise UserError("This protocol is limited to 96 sequenced colonies, please "
                        "submit additional runs if needed.")
    cols = int(math.ceil(total_cols / float(8)))
    columns = [{"column": i, "volume": "150:microliter"} for i in range(0, cols)]
    protocol.dispense(growth_plate, params['growth_media'], columns)

    growth_wells = growth_plate.all_wells(columnwise=True)

    for i, c in enumerate(constructs):
        protocol.autopick(agar_wells[i],
                          growth_wells[i*num_colonies:i*num_colonies+num_colonies],
                          #min_count=0,
                          dataref="%s" % agar_wells[i].name)
        for j in range(num_colonies):
            growth_wells[i*num_colonies+j].set_name('%s_growth_%s' % (agar_wells[i].name, j+1))
            if 'mutant_objs' in params.keys():
                mut = next(m for m in params['mutant_objs'] if m.name == agar_wells[i].name)
                mut.add_growth_wells(growth_wells[i*num_colonies+j])

    protocol.cover(growth_plate, lid="low_evaporation")
    protocol.incubate(growth_plate, "warm_37", "16:hour", shaking=True, co2=0)

    # return growth plate for uncovering in full script
    return growth_plate

if __name__ == '__main__':
    from autoprotocol.harness import run
    run(transform, 'Transform')
