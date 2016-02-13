from autoprotocol.instruction import Instruction

'''
    :copyright: 2015 by The Autoprotocol Development Team, see AUTHORS
        for more details.
    :license: BSD, see LICENSE for more details

'''


class Miniprep(Instruction):
    """
    Perform mini-prep on indicated aliquot(s).

    Parameters
    ----------
    group_opts : list
        List of {from: well1, to: well2} instructions

    """
    def __init__(self, group_opts):
        super(Miniprep, self).__init__({
            "op": "miniprep",
            "groups": group_opts
        })
