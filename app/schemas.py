from typing import List, Optional

from pydantic import BaseModel


class PredictionInput(BaseModel):
    smiles: List[str] = []
    inchis: List[str] = []


class Prediction(BaseModel):
    atom_id: int
    min: float
    mean: float
    max: float


class PredictionResponse(BaseModel):
    input: str
    standard_smiles: str
    molblock: str
    predictions: List[Prediction]
