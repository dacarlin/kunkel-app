from autoprotocol.container import WellGroup
from autoprotocol import UserError
from autoprotocol_lib import transcriptic_protocol
from modules.utils import return_dispense_media, printdatetime
import math

def plasmidprep(protocol, params):
    names = [item['name'] for item in params['samples']]
    num_samples = len(params['samples'])
    if len(names) != len(set(names)):
        raise UserError("All plasmidprep names have to be unique.")
    if num_samples > 96:
        raise UserError("A maximum of 96 preps per run is allowed.")
    if params['type'] == "Miniprep":
        prep = protocol.miniprep
        dest_plate = protocol.ref('miniprep_%s' % printdatetime(),
                                  cont_type="96-pcr", storage="cold_20")
        transfer_vol = "35:microliter"
        duration = "16:hour"
    elif params['type'] == 'Maxiprep':
        prep = protocol.maxiprep
        transfer_vol = "300:microliter"
        dest_plate = protocol.ref('maxiprep_%s' % printdatetime(),
                                  cont_type="96-deep", storage="cold_20")
        duration = "16:hour"

    dests = dest_plate.wells_from(0, num_samples, columnwise=True)
    cols = int(math.ceil(float(num_samples)/8.0))
    dispense_cols = [{"column": col, "volume": "1800:microliter"}
                     for col in range(0, cols)]

    # Aliquot, grow, prep
    growth_plate = protocol.ref("plasmid_prep_growth_plate_%s" % printdatetime(),
                                 cont_type="96-deep", discard=True)
    protocol.dispense(growth_plate, params['media'], dispense_cols)
    src_samples = WellGroup([item['sample'] for item in params['samples']])
    protocol.transfer(src_samples,
                      growth_plate.wells_from(0, num_samples, columnwise=True),
                      "10:microliter")

    protocol.cover(params["growth_plate"], lid="low_evaporation")
    protocol.cover(growth_plate, lid="standard")
    protocol.incubate(growth_plate, "warm_37", duration, shaking=True, co2=0)

    prep(growth_plate.wells_from(0, num_samples, columnwise=True),
         dests)

    if params['type'] == 'Maxiprep':
        new_dests = []
        for i in range(0, num_samples):
            new_dests.append(protocol.ref('maxiprep_%d_%s' % (i+1, printdatetime()),
                                      cont_type="micro-1.5",
                                      storage="cold_20").well(0))
            protocol.transfer(dests[i], new_dests[i], transfer_vol)
        dests = new_dests

    for name, well in zip(names, dests):
        well.set_name("%s_%s" % (name, printdatetime(time=False)))

if __name__ == '__main__':
    from autoprotocol.harness import run
    run(plasmidprep, "PlasmidPrep")
