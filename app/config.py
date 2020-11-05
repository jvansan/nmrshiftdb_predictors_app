import os
from pathlib import Path

CWD = Path(os.path.realpath(__file__)).parent
LIB_DIR = CWD.joinpath("lib")
TEMP_DIR = CWD.joinpath("temp")
PRED_C_JAR = "predictorc.jar:cdk-2.3.jar"
PRED_H_JAR = "predictorh.jar:cdk-2.3.jar"
PREFIX = "/api"
