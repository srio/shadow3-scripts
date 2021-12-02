
from syned.util.json_tools import load_from_json_file

from syned.storage_ring.electron_beam import ElectronBeam
from syned.storage_ring.magnetic_structures.undulator import Undulator
from syned.beamline.optical_elements.ideal_elements.screen import Screen
from syned.beamline.optical_elements.ideal_elements.lens import IdealLens
from syned.beamline.optical_elements.absorbers.filter import Filter
from syned.beamline.optical_elements.absorbers.slit import Slit
from syned.beamline.optical_elements.absorbers.beam_stopper import BeamStopper
from syned.beamline.optical_elements.mirrors.mirror import Mirror
from syned.beamline.optical_elements.crystals.crystal import Crystal
from syned.beamline.optical_elements.gratings.grating import Grating

from syned.beamline.shape import SurfaceShape, Conic, Ellipsoid, Plane
from syned.beamline.shape import Rectangle
from syned.storage_ring.light_source import LightSource

from syned.beamline.beamline import Beamline
from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates



if True:
    #
    beamline = load_from_json_file("/users/srio/Oasys/id03.json")

    # print(beamline.info())

    # print(tmp.get_light_source().info())
    element = ["CPMU18"]
    indices = [0]
    dist_to_source = [0.0]
    type1 = ['source']
    h = [0.0]
    v = [0.0]
    thickness = [None]
    formula = [None]
    density = [None]

    dist_cumulated = 0.0
    for i,element_i in enumerate(beamline.get_beamline_elements()):
        oe_i = element_i.get_optical_element()

        element.append(oe_i.get_name())
        indices.append(i)

        coor = element_i.get_coordinates()
        dist_cumulated += coor.p()
        dist_to_source.append(dist_cumulated)

        if isinstance(oe_i, Filter):
            type1.append('window')
        elif isinstance(oe_i, Slit):
            type1.append('slit')
        else:
            type1.append('unknown')


        shape = oe_i.get_boundary_shape()
        if isinstance(shape, Rectangle):
            x_left, x_right, y_bottom, y_top = shape.get_boundaries()
            h.append(x_right - x_left)
            v.append(y_top - y_bottom)
        else:
            h.append(None)
            v.append(None)

        if isinstance(oe_i, Filter):
            thickness.append(oe_i.get_thickness() * 1e3)
        else:
            thickness.append(None)

        if isinstance(oe_i, Filter):
            formula.append(oe_i.get_material())
        else:
            formula.append(None)

        if isinstance(oe_i, Filter):
            if oe_i.get_material() == "C":
                density.append(3.52)
            else:
                density.append(None)
        else:
            density.append(None)





    print("element: ",element)
    print("dist_to_source", dist_to_source)
    print("type: ", type1)
    print("h: ",h)
    print("v: ",v)
    print("thickness: ",thickness)
    print("formula: ", formula)
    print("density: ", density)

    #
    # TODO
    #
    # implement name
    # filter does not have size
    # implement in Oasys automatic propagation
    # implement scripr widget

