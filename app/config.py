import os
from pathlib import Path

CWD = Path(os.path.realpath(__file__)).parent
LIB_DIR = CWD.joinpath("lib")
TEMP_DIR = CWD.joinpath("temp")
PRED_C_JAR = "predictorc.jar:cdk-interfaces.jar:cdk-io.jar:cdk-core.jar:cdk-data.jar:cdk-standard.jar:cdk-sdg.jar:cdk-valencycheck.jar:cdk-atomtype.jar:cdk-legacy.jar:cdk-smiles.jar:cdk-ctab.jar:vecmath1.2-1.14.jar:jgrapht.jar:guava-17.0.jar:beam-core-0.9.2.jar"
PRED_H_JAR = "predictorh.jar:cdk-interfaces.jar:cdk-io.jar:cdk-core.jar:cdk-data.jar:cdk-standard.jar:cdk-sdg.jar:cdk-valencycheck.jar:cdk-atomtype.jar:cdk-legacy.jar:cdk-smiles.jar:cdk-ctab.jar:vecmath1.2-1.14.jar:jgrapht.jar:guava-17.0.jar:beam-core-0.9.2.jar"
PREFIX = "/api"
