from autoprotocol.container import Container, WellGroup
from autoprotocol.protocol import Ref
from autoprotocol.unit import Unit
import datetime


def find_named_well(name, wells):
    for w in wells.wells:
        if w.name == name:
            return w

def printdatetime(time=True):
    printdate = datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
    if not time:
        printdate = datetime.datetime.now().strftime('%Y-%m-%d')
    return printdate

def provision_to_tube(protocol, name, tube, resource_id, volume, discard=True, storage=None):
    '''
        provision_to_tube allows us to provision into a tube/well that can then be used
        in transfers on the workcell.
    '''
    assert isinstance(volume, (Unit, int, float)), "Volume must be of type int, float or Unit."
    if isinstance(volume, Unit):
        volume = volume.value
    if storage:
        dest = protocol.ref(name, None, tube, storage=storage).well(0)
    else:
        dest = protocol.ref(name, None, tube, discard=discard).well(0)
    protocol.provision(resource_id, dest, "%s:microliter" % volume)
    return(dest)


def ref_kit_container(protocol, name, container, kit_id, discard=True, store=None):
    '''
        Still in use to allow booking of agar plates on the fly
    '''
    kit_item = Container(id, protocol.container_type(container), name=name)
    if store:
        protocol.refs[name] = Ref(name, {"reserve": kit_id, "store": {"where": store}}, kit_item)
    else:
        protocol.refs[name] = Ref(name, {"reserve": kit_id, "discard": discard}, kit_item)
    return(kit_item)


def make_list(my_str, integer=False):
    '''
        Sometimes you need a list of a type that is not supported. This takes a string and
        comma-seperates it, returning a list of strings or integers.
    '''
    assert isinstance(my_str, str), "Input needs to be of type string."
    if integer:
        my_str = [int(x.strip()) for x in my_str.split(",")]
    else:
        my_str = [x.strip() for x in my_str.split(",")]
    return my_str


def thermocycle_ramp(start_temp, end_temp, total_duration, step_duration):
    '''
        Create a ramp instruction for the thermocyler. Used in annealing protocols.
    '''
    assert Unit.fromstring(total_duration).unit == Unit.fromstring(step_duration).unit, ("Thermocycle_ramp durations"
                                                                                         " must be specified using the"
                                                                                         " same unit of time.")
    thermocycle_steps = []
    start_temp = Unit.fromstring(start_temp).value
    num_steps = int(Unit.fromstring(total_duration).value // Unit.fromstring(step_duration).value)
    step_size = (Unit.fromstring(end_temp).value - start_temp) // num_steps
    for i in xrange(0, num_steps):
        thermocycle_steps.append({
            "temperature": "%d:celsius" % (start_temp + i * step_size), "duration": step_duration
            })
    return thermocycle_steps


def return_agar_plates(wells):
    '''
        Dicts of all plates available that can be purchased.
    '''
    if wells == 6:
        plates = {"lb-broth-50ug-ml-kan": "ki17rs7j799zc2",
                  "lb-broth-100ug-ml-amp": "ki17sbb845ssx9",
                  "lb-broth-100ug-ml-specto": "ki17sbb9r7jf98",
                  "lb-broth-100ug-ml-cm": "ki17urn3gg8tmj",
                  "noAB": "ki17reefwqq3sq"}
    elif wells == 1:
        plates = {"lb-broth-50ug-ml-kan": "ki17t8j7kkzc4g",
                  "lb-broth-100ug-ml-amp": "ki17t8jcebshtr",
                  "lb-broth-100ug-ml-specto": "ki17t8jaa96pw3",
                  "lb-broth-100ug-ml-cm": "ki17urn592xejq",
                  "noAB": "ki17t8jejbea4z"}
    else:
        raise ValueError("Wells has to be an integer, either 1 or 6")
    return(plates)


def det_new_group(i, base=0):
    '''
        Helper to determine if new_group should be added. Returns true when i matches the base, which defaults to 0.
    '''
    assert isinstance(i, int), "Needs an integer."
    assert isinstance(base, int), "Base has to be an integer"
    if i == base:
        new_group = True
    else:
        new_group = False
    return new_group


def return_dispense_media():
    '''
        Dict of media for reagent dispenser
    '''
    media = {"50_ug/ml_Kanamycin": "lb-broth-50ug-ml-kan",
             "100_ug/ml_Ampicillin": "lb-broth-100ug-ml-amp",
             "100_ug/mL_Spectinomycin": "lb-broth-100ug-ml-specto",
             "30_ug/ml_Kanamycin": "lb-broth-30ug-ml-kan",
             "50_ug/ml_Kanamycin_25_ug/ml_Chloramphenicol": "lb-broth-50ug-ml-kan-25ug-ml-cm",
             "15_ug/ml_Tetracycline": "lb-broth-15ug-ml-tet",
             "25_ug/ml_Chloramphenicol": "lb-broth-25ug-ml-cm",
             "LB_broth": "lb-broth-noAB"}
    return(media)
