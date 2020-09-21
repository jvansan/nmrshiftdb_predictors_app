from typing import List, Optional
from enum import Enum

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


class CanonMethods(str, Enum):
    smiles = "smiles"
    algorithm = "algorithm"


class CanonInput(BaseModel):
    smiles: str
    add_hs: bool
    reverse: bool
    method: CanonMethods


class CanonOutput(BaseModel):
    input: CanonInput
    result: str
