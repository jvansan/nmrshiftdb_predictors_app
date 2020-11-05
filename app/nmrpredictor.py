import subprocess
from pathlib import Path
from enum import Enum
from typing import Optional
from app import config


class Predictors(Enum):
    C = 1
    H = 2

class Solvent(str, Enum):
    cdcl3 = "Chloroform-D1 (CDCl3)"
    cd3od = "Methanol-D4 (CD3OD)"
    c2d6s0 = "Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)"

def run_predictorc(path: Path, solvent: Optional[Solvent]):
    if not path.exists():
        raise FileNotFoundError(f"Missing file {path} - cannot run predictor")
    res = call_predictor(path, Predictors.C, solvent=solvent)
    return parse_predictor_output(res)


def run_predictorh(path: Path, solvent: Optional[Solvent]):
    if not path.exists():
        raise FileNotFoundError(f"Missing file {path} - cannot run predictor")
    res = call_predictor(path, Predictors.H, solvent=solvent)
    return parse_predictor_output(res)


def call_predictor(path: Path, predictor: Predictors, solvent: Optional[Solvent]):
    if predictor == Predictors.C:
        jar = config.PRED_C_JAR
    elif predictor == Predictors.H:
        jar = config.PRED_H_JAR
    else:
        raise ValueError("Predictor does not exist")
    if solvent:
        cmd = ["java", "-cp", jar, "Test", str(path.absolute()), f'\"{solvent.value}\"' ]
    else:
        cmd = ["java", "-cp", jar, "Test", str(path.absolute())]
    return subprocess.run(" ".join(cmd), cwd=config.LIB_DIR, capture_output=True, shell=True).stdout


def parse_predictor_output(res):
    res_list = res.decode().split("\n")
    data = []
    for line in res_list[2:]:
        x = line.split()
        if not x:
            continue
        data.append(
            {
                "atom_id": int(x[0].replace(":", "")),
                "min": float(x[1]),
                "mean": float(x[2]),
                "max": float(x[3]),
            }
        )
    return data
