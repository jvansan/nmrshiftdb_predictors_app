# nmrshiftdb_predictors_app

Web API wrapper for [nmrshiftdb2 predictors](https://sourceforge.net/p/nmrshiftdb2/wiki/PredictorJars/).
Please see https://nmrshiftdb.nmr.uni-koeln.de for details on nmrshiftdb.

This app will parse a list of InChI and SMILES strings and return canonicalized atom numbering with predictions.

The <sup>1</sup>H predictor runs much faster than the <sup>13</sup>C predictor.

nmrshiftdb and [CDK](https://cdk.github.io/) JARs are distributed in the `app/lib` directory.

nmrshiftdb is licensed under a GNU AGPL license.

CDK is licensed under a GNU LGLP 2.1 license.

#### Setup

Requires JDK, and Anaconda Python 3.

To install dependencies (only for running locally):

```
conda install -n nmr_api -c conda-forge fastapi uvicorn rdkit
```

#### To run locally

`conda activate nmr_api && uvicorn app.main:app --reload --loop="asyncio"`

#### To run in Docker

First build the docker image:

`docker build -t nmr_api:latest .`

Then run the docker image:

`docker run --name nmr_api -t -d -p 8000:80 nmr_api`

And then you can go to `http://localhost:8000/api/docs` and play with the API docs.

#### Usage in Python

Here is an example of how to use it from Python:

Possible solvents are "Chloroform-D1 (CDCl3)", "Methanol-D4 (CD3OD)", and "Dimethylsulphoxide-D6 (DMSO-D6, C2D6SO)".
Note you will have the url-encode these strings using the [urlib.parse](https://docs.python.org/3/library/urllib.parse.html#url-quoting) library.

Excluded the solvent to use the "Unreported" option:
> The option "Unreported" will use shifts from any solvent.

```
import requests
from urllib.parse import quote
from rdkit import Chem
smiles = ["C", "CC", "CCC", "c1ccccc1"]
solvent = "Chloroform-D1 (CDCl3)"
r = requests.post("http://localhost:8000/api/predict/proton?solvent="+quote(solvent), json={"smiles": smiles})
proton_preds = r.json()
r1 = requests.post("http://localhost:8000/api/predict/carbon?solvent="+quote(solvent), json={"smiles": smiles})
carbon_preds = r1.json()

benzene_prot = proton_preds[3]
benzene_carb = carbon_preds[3]

# you can write molfiles from response for visualization in ChemDraw
with open("benzene.mol", "w") as f:
    f.write(benzene_prot["molblock"])

# you can load the responses into RDKit and write images
from rdkit.Chem import Draw

# Don't sanitize in order to keep Hs in parsing, then sanitize to cleanup errors when writing image
m = Chem.MolFromSmiles(benzene_prot['standard_smiles'], sanitize=False)
Chem.SantizeMol(m)
Draw.MolToImageFile(m, "benzene.png", includeAtomNumbers=True)
```
