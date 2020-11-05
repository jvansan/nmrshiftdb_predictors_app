from typing import List, Optional

from fastapi import FastAPI
from starlette.responses import RedirectResponse

from app import schemas, chem
from app.chem import process_inchi, process_smiles
from app.config import PREFIX
from app.nmrpredictor import run_predictorc, run_predictorh, Solvent

app = FastAPI(
    title="NMRshiftDB2 predictor app",
    redoc_url=None,
    docs_url=f"{PREFIX}/docs",
    openapi_url=f"{PREFIX}/openapi.json",
)


@app.get(f"{PREFIX}/status", tags=["Root"], include_in_schema=False)
async def status():
    return {"status": "good"}


@app.get("/", tags=["Root"], include_in_schema=False)
@app.get(f"{PREFIX}/", tags=["Root"], include_in_schema=False)
async def doc_redirect():
    return RedirectResponse(url=f"{PREFIX}/docs")


@app.post(f"{PREFIX}/predict/proton", response_model=List[schemas.PredictionResponse])
def predict_proton(input_: schemas.PredictionInput, solvent: Optional[Solvent] = None):
    """Predicts 1H NMR peaks for input SMILES and INCHI"""
    res = []
    solv = "Unreported"
    if solvent is not None:
        solv = solvent.value 
    for smi in input_.smiles:
        canon_smi, canon_mb, temp_path = process_smiles(smi)
        pred = run_predictorh(temp_path, solvent=solvent)
        res.append(
            schemas.PredictionResponse(
                input=smi,
                solvent=solv,
                standard_smiles=canon_smi,
                molblock=canon_mb,
                predictions=pred,
            )
        )
    for inchi in input_.inchis:
        canon_smi, canon_mb, temp_path = process_inchi(inchi)
        pred = run_predictorh(temp_path)
        res.append(
            schemas.PredictionResponse(
                input=smi,
                solvent=solv,
                standard_smiles=canon_smi,
                molblock=canon_mb,
                predictions=pred,
            )
        )
    return res


@app.post(f"{PREFIX}/predict/carbon", response_model=List[schemas.PredictionResponse])
def predict_carbon(input_: schemas.PredictionInput, solvent: Optional[Solvent] = None):
    """Predicts 13C NMR peaks for input SMILES and INCHI"""
    res = []
    solv = "Unreported"
    if solvent is not None:
        solv = solvent.value 
    for smi in input_.smiles:
        canon_smi, canon_mb, temp_path = process_smiles(smi)
        pred = run_predictorc(temp_path, solvent=solvent)
        res.append(
            schemas.PredictionResponse(
                input=smi,
                standard_smiles=canon_smi,
                molblock=canon_mb,
                solvent=solv,
                predictions=pred,
            )
        )
    for inchi in input_.inchis:
        canon_smi, canon_mb, temp_path = process_inchi(inchi)
        pred = run_predictorc(temp_path)
        res.append(
            schemas.PredictionResponse(
                input=smi,
                standard_smiles=canon_smi,
                solvent=solv,
                molblock=canon_mb,
                predictions=pred,
            )
        )
    return res


@app.post(f"{PREFIX}/utils/canonicalize", response_model=schemas.CanonOutput)
def canonicalize_smiles(
    smiles: str,
    method: schemas.CanonMethods = schemas.CanonMethods.smiles,
    reverse: Optional[bool] = False,
    add_hs: Optional[bool] = True,
):
    """Perform canonicalization as defined by the input.
    Reverse only valid for algorithm implementation."""
    m = chem.smi_to_mol(smiles)
    if method == schemas.CanonMethods.algorithm:
        result = chem.canonicalize_atom_order(m, reverse=reverse, add_hs=add_hs)
    elif method == schemas.CanonMethods.smiles:
        result = chem.rdkit_atom_order(m, add_hs=add_hs)
    else:
        raise TypeError("Invalid input method")
    return {
        "input": {
            "smiles": smiles,
            "method": method,
            "reverse": reverse,
            "add_hs": add_hs,
        },
        "result": chem.mol_to_smi(result),
    }
