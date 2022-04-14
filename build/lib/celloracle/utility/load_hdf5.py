# -*- coding: utf-8 -*-



from ..network_analysis import load_links
from ..network import load_net
from ..trajectory.oracle_core import load_oracle
from ..motif_analysis.tfinfo_core import load_TFinfo
from ..applications.differentiation_flow import load_gradient

def load_hdf5(file_path, object_class_name=None):
    """
    Load an object of celloracle's custom class that was saved as hdf5.

    Args:
        file_path (str): file_path.
        object_class_name (str): Types of object.
            If it is None, object class will be identified from the extension of file_name.
            Default is None.
    """

    # identify object class type
    if object_class_name is None:
        object_class_name = file_path.split(".")[-1]
    object_class_name = object_class_name.capitalize()

    # load object
    if object_class_name == "Links":
        obj = load_links(file_path=file_path)

    elif object_class_name == "Net":
        obj = load_net(file_path=file_path)

    elif object_class_name == "Oracle":
        obj = load_oracle(file_path=file_path)

    elif object_class_name in ["Tfinfo", "tfinfo"]:
        obj = load_TFinfo(file_path=file_path)

    elif object_class_name == "Gradient":
        obj = load_gradient(file_path=file_path)

    else:
        print(f"object_class_name: {object_class_name} was not in the loading option")
        raise ValueError("File type identification failed. Please enter 'object_class_name' manually.")

    return obj
