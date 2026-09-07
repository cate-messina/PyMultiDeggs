import pandas as pd
import pytest
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def rnaseq():
    return pd.read_csv(DATA_DIR / "synthetic_rnaseqData.csv", index_col=0)


@pytest.fixture
def proteomics():
    return pd.read_csv(DATA_DIR / "synthetic_proteomicData.csv", index_col=0)


@pytest.fixture
def olink():
    return pd.read_csv(DATA_DIR / "synthetic_OlinkData.csv", index_col=0)


@pytest.fixture
def metadata():
    return pd.read_csv(DATA_DIR / "synthetic_metadata.csv", index_col=0)


@pytest.fixture
def assayData_multiomic(rnaseq, proteomics, olink):
    return {
        "RNAseq": rnaseq,
        "Proteomics": proteomics,
        "Olink": olink,
    }