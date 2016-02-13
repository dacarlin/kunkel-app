from autoprotocol.protocol import Protocol
from autoprotocol.container import WellGroup
from transcriptic_instruction import *
import sys
if sys.version_info.major == 3:
    basestring = str

'''
    :copyright: 2016 by The Autoprotocol Development Team, see AUTHORS
        for more details.
    :license: MIT, see LICENSE for more details

'''

"""
    NOTE : Import order matters when using this extension
    Do `from autoprotocol import Protocol` BEFORE
    `import transcriptic_protocol`
"""


def miniprep(self, source_well, dest_well):
    """
    Perform mini-prep on source wells and resulting DNA aliquot will be placed
    in destination wells.

    Supports a one-to-one and many-to-one mapping of source to destination
    wells.

    Example Usage:

    .. code-block:: python

        p = Protocol()
        my_cell_plate = p.ref("my_cell_plate", None, "96-flat", discard = True)
        my_dna_plate = p.ref("my_dna_plate", None, "96-flat", storage = "cold_4")

        p.miniprep(my_cell_plate.wells_from(0,3), my_dna_plate.well(0))

    Autoprotocol Output:

    .. code-block:: json

        "instructions": [
            {
              "op": "miniprep",
              "groups": [{
                "from": "my_cell_plate/0",
                "to": "my_dna_plate/0"
              },{
                "from": "my_cell_plate/1",
                "to": "my_dna_plate/0"
              },{
                "from": "my_cell_plate/2",
                "to": "my_dna_plate/0"
              }]
            }
          ]

    Parameters
    ----------
    source_well : Well, WellGroup
      Well(s) containing bacteria to be mini-prepped
    dest_well : Well, WellGroup
      Well(s) indicating final storage location(s) of resulting DNA

    """
    source = WellGroup(source_well)
    dest = WellGroup(dest_well)
    len_source = len(source.wells)
    len_dest = len(dest.wells)

    if len_source != len_dest and (len_dest != 1):
        raise ValueError("Number of source wells need to match number of "
                         "destination wells unless only one destination well "
                         "is specified.")

    if len_source > 1 and len_dest == 1:
        dest = WellGroup(dest.wells * len_source)

    group_opts = []
    for s_well, d_well in zip(source.wells, dest.wells):
        group_dict = {}
        group_dict["from"] = s_well
        group_dict["to"] = d_well
        group_opts.append(group_dict)

    self.instructions.append(Miniprep(group_opts))

old_as_dict = Protocol.as_dict
def as_dict(self):
    res = old_as_dict(self)
    if hasattr(self, "price"):
        res['price'] = self.price
    return res


'''Bind new methods to existing protocol'''
Protocol.miniprep = miniprep
Protocol.as_dict = as_dict
