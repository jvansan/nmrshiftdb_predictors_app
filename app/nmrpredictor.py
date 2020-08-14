import subprocess
from pathlib import Path
from enum import Enum
from app import config


class Predictors(Enum):
    C = 1
    H = 2


def run_predictorc(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file {path} - cannot run predictor")
    res = call_predictor(path, Predictors.C)
    return parse_predictor_output(res)


def run_predictorh(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file {path} - cannot run predictor")
    res = call_predictor(path, Predictors.H)
    return parse_predictor_output(res)


def call_predictor(path: Path, predictor: Predictors):
    if predictor == Predictors.C:
        jar = config.PRED_C_JAR
    elif predictor == Predictors.H:
        jar = config.PRED_H_JAR
    else:
        raise ValueError("Predictor does not exist")
    return subprocess.check_output(
        ["java", "-cp", jar, "Test", str(path.absolute())], cwd=config.LIB_DIR,
    )


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
