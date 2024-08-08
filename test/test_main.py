import pytest
from fastapi import HTTPException
from src.main import substructure_search


# Check positive scenarios
@pytest.mark.parametrize(
    "substructure,expected",
    [("CCO", [{'molecule_id': 2, 'smiles_structure': 'CCO'}, {'molecule_id': 3, 'smiles_structure': 'CC(=O)O'}, {'molecule_id': 4, 'smiles_structure': 'CC(=O)Oc1ccccc1C(=O)O'}]),
     ("c1ccccc1", [{'molecule_id': 1, 'smiles_structure': 'Cc1ccccc1'}, {'molecule_id': 4, 'smiles_structure': 'CC(=O)Oc1ccccc1C(=O)O'}]),
     ("CC(=O)Oc1ccccc1C(=O)O",[{'molecule_id': 4, 'smiles_structure': 'CC(=O)Oc1ccccc1C(=O)O'}])])
def test_substructure_search_pass(substructure, expected):
    """Test positive scenarios where the substructure should be found."""
    result = substructure_search(substructure)
    assert result == expected

# Check negative scenarios - substructures for which there are no matching structures
@pytest.mark.xfail
@pytest.mark.parametrize(
    "substructure,expected",
    [("O=C=O", [{'molecule_id': 1, 'smiles_structure': 'Cc1ccccc1'}, {'molecule_id': 4, 'smiles_structure': 'CC(=O)Oc1ccccc1C(=O)O'}]),
     ("FC(Br)(Cl)F", [{'molecule_id': 4, 'smiles_structure': 'CC(=O)Oc1ccccc1C(=O)O'}])])
def test_substructure_search_fail(substructure, expected):
    """Test negative scenarios where no matches are expected for the substructure."""
    result = substructure_search(substructure)
    assert result == expected

# Check scenarios where invalid SMILES structures should raise an HTTPException
@pytest.mark.parametrize(
    "not_smiles_structure",
    ["invalid", "1234", "!@#$"])
def test_substructure_search_fail_error_expected(not_smiles_structure):
    """Test scenarios where an invalid SMILES structure should raise an HTTPException."""
    with pytest.raises(HTTPException):
        substructure_search(not_smiles_structure)







