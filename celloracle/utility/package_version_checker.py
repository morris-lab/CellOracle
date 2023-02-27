
import os
import pandas as pd
#from warnings import warn

try:
    from importlib.metadata import version
except:
    from importlib_metadata import version

from ..utility import __path__ as PATH


def _load_requirements():
    path = os.path.abspath(PATH[0])
    REQUIREMENTS_PATH = os.path.join(path, "requirements.txt")

    with open(REQUIREMENTS_PATH) as f:
        requirements = f.read().splitlines()

    return requirements


def _get_package_names(requirements):
    names = []
    for i in requirements:
        if ">" in i:
            i = i.split(">")[0]
        elif "=" in i:
            i = i.split("=")[0]
        names.append(i)

    return names

def _get_installed_version(requirements):
    names = _get_package_names(requirements)
    installed_versions = []
    for i in names:
        try:
            v = version(i)
        except:
            v = "not_found"
        installed_versions.append(v)

    return installed_versions

def _get_required_version(requirements):
    versions = []
    for i in requirements:
        if ">=" in i:
            i = i.split(">=")[-1]
        elif ">" in i:
            i = i.split(">")[-1]
        elif "=" in i:
            i = i.split("=")[-1]
        else:
            i = "auto"
        versions.append(i)

    return versions

def _is_version_OK(que, ref):

    try:
        answer = True
        if ref == "auto":
            if que == "not_found":
                answer = False
            else:
                pass
        else:

            que_ = que.split(".")
            ref_ = ref.split(".")


            for q, r in zip(que_, ref_):
                #print(q, r)
                if int(q) < int(r):
                    answer = False
                    break
                elif int(q) > int(r):
                    answer = True
                    break
                else:
                    pass
        return answer
    except:
        return "unknown"

def check_python_requirements(return_detail=True, print_warning=True):
    """
    Check installation status and requirements of dependant libraries.
    """

    try:
        REQUIREMENTS = _load_requirements()
    except:
        if print_warning:
            print("Could not check requirements.")
        return None

    status = pd.DataFrame({"package_name": _get_package_names(REQUIREMENTS),
                           "installed_version": _get_installed_version(REQUIREMENTS),
                            "required_version": _get_required_version(REQUIREMENTS)})
    status["requirement_satisfied"] = \
    [_is_version_OK(ref, que) for ref, que in zip(status["installed_version"],
                                                 status["required_version"])]

    if print_warning:
        n_not_ok = (status["requirement_satisfied"] == False).sum()
        if n_not_ok == 1:
            print(f"{n_not_ok} package does not meet CellOracle requirement.")
        elif n_not_ok >= 2:
            print(f"{n_not_ok} packages do not meet CellOracle requirement.")

        for i, (name, installed, required, isok) in status.iterrows():
            if isok is False:
                print(f" Your {name} version is {installed}. Please install {REQUIREMENTS[i]}")
    if return_detail:
        return status
