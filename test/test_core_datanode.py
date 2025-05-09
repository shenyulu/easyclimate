"""
pytest for core/datanode.py
"""

import pytest
from easyclimate.core.datanode import DataNode, open_datanode
import xarray as xr
import numpy as np
import json
import shutil


@pytest.fixture
def sample_node():
    """Fixture to create a sample DataNode hierarchy for testing"""
    node = DataNode(name="root")
    node.a = 42
    node.b = "test"
    node.subnode = DataNode(name="sub")
    node.subnode.x = 3.14
    node.subnode.y = xr.DataArray(np.random.rand(3, 3), name="test_array")
    return node


def test_initialization():
    """Test DataNode initialization"""
    node = DataNode(name="test")
    assert node.name == "test"
    assert node._attributes == {}


def test_attribute_access(sample_node):
    """Test attribute access and auto-creation"""
    assert sample_node.a == 42
    assert sample_node.b == "test"
    assert isinstance(sample_node.subnode, DataNode)
    assert sample_node.subnode.x == 3.14

    # Test auto-creation of new nodes
    assert isinstance(sample_node.new_node, DataNode)
    assert sample_node.new_node.name == "new_node"


def test_item_access(sample_node):
    """Test dictionary-style access"""
    assert sample_node["a"] == 42
    assert sample_node["b"] == "test"
    assert sample_node["subnode"]["x"] == 3.14

    # Test path-style access
    assert sample_node["subnode/x"] == 3.14


def test_item_assignment(sample_node):
    """Test dictionary-style assignment"""
    sample_node["c"] = 100
    assert sample_node.c == 100

    # Test path-style assignment
    sample_node["subnode/new"] = "value"
    assert sample_node.subnode.new == "value"

    # Test nested path creation
    sample_node["path/to/new/node"] = 99
    assert sample_node.path.to.new.node == 99


def test_repr(sample_node):
    """Test string representation"""
    repr_str = repr(sample_node)
    assert "root" in repr_str
    assert "a: 42" in repr_str
    assert "sub" in repr_str


def test_html_representation(sample_node):
    """Test HTML representation generation"""
    html = sample_node._repr_html_()
    assert isinstance(html, str)
    assert "root" in html
    assert "sub" in html
    assert "xarray.DataArray" in html


def test_xarray_integration(sample_node):
    """Test xarray object handling"""
    assert isinstance(sample_node.subnode.y, xr.DataArray)
    assert sample_node.subnode.y.shape == (3, 3)

    # Test HTML representation contains xarray's HTML
    html = sample_node._repr_html_()
    assert "dim_0" in html  # xarray default dimension name


def test_zarr_metadata(tmp_path, sample_node):
    """Test zarr metadata structure"""
    save_path = tmp_path / "test_node"
    sample_node.to_zarr(save_path)

    # Check root metadata
    with open(save_path / "metadata.json") as f:
        metadata = json.load(f)

    assert metadata["name"] == "root"
    assert metadata["attributes"]["a"] == 42
    assert isinstance(metadata["attributes"]["subnode"], dict)
    assert metadata["attributes"]["subnode"]["__type__"] == "DataNode"

    # Check subnode metadata
    with open(save_path / "subnode" / "metadata.json") as f:
        sub_metadata = json.load(f)

    assert sub_metadata["name"] == "sub"
    assert sub_metadata["attributes"]["x"] == 3.14
    assert sub_metadata["attributes"]["y"]["__type__"] == "xarray"


def test_format_value():
    """Test value formatting for display"""
    node = DataNode(name="test")

    # Test list formatting
    long_list = list(range(10))
    formatted = node._format_value(long_list)
    assert "[0, 1, ..., 9]" in formatted
    assert "length: 10" in formatted

    # Test dict formatting
    long_dict = {i: i * 2 for i in range(10)}
    formatted = node._format_value(long_dict)
    assert "0: ..." in formatted
    assert "length: 10" in formatted

    # Test simple value
    assert node._format_value(42) == "42"


def test_open_datanode_basic(tmp_path, sample_node):
    """Test basic loading functionality of open_datanode"""
    save_path = tmp_path / "test_node"
    sample_node.to_zarr(save_path)

    # Load using open_datanode
    loaded_node = open_datanode(save_path)

    assert loaded_node.name == "root"
    assert loaded_node.a == 42
    assert loaded_node.b == "test"
    assert isinstance(loaded_node.subnode, DataNode)
    assert loaded_node.subnode.x == 3.14


def test_open_datanode_nested_structure(tmp_path, sample_node):
    """Test loading nested node structure"""
    save_path = tmp_path / "nested_node"
    sample_node.to_zarr(save_path)

    loaded_node = open_datanode(save_path)

    assert isinstance(loaded_node.subnode, DataNode)
    assert loaded_node.subnode.name == "sub"
    assert loaded_node.subnode.x == 3.14


def test_open_datanode_path_types(tmp_path, sample_node):
    """Test that open_datanode works with both str and Path inputs"""
    save_path = tmp_path / "path_type_node"
    sample_node.to_zarr(save_path)

    # Test with Path object
    loaded_path = open_datanode(save_path)
    # Test with string path
    loaded_str = open_datanode(str(save_path))

    assert loaded_path.name == loaded_str.name
    assert loaded_path.subnode.x == loaded_str.subnode.x
    assert loaded_path.subnode.y.equals(loaded_str.subnode.y)


def test_open_datanode_invalid_path():
    """Test error handling for invalid paths"""
    with pytest.raises(FileNotFoundError):
        open_datanode("nonexistent/path")


def test_open_datanode_missing_metadata(tmp_path):
    """Test error handling when metadata is missing"""
    invalid_path = tmp_path / "invalid_node"
    invalid_path.mkdir()

    with pytest.raises(FileNotFoundError):
        open_datanode(invalid_path)


def test_open_datanode_corrupted_metadata(tmp_path, sample_node):
    """Test error handling for corrupted metadata"""
    save_path = tmp_path / "corrupted_node"
    sample_node.to_zarr(save_path)

    # Corrupt the metadata file
    with open(save_path / "metadata.json", "w") as f:
        f.write("{invalid json")

    with pytest.raises(json.JSONDecodeError):
        open_datanode(save_path)


def test_open_datanode_complex_hierarchy(tmp_path):
    """Test loading a complex node hierarchy"""
    # Create complex structure
    node = DataNode(name="complex")
    node.num = 123
    node.text = "hello"
    node.data = xr.DataArray(np.random.rand(2, 2), name="data")
    node.group = DataNode(name="group")
    node.group.a = 1
    node.group.b = 2
    node.group.sub = DataNode(name="sub")
    node.group.sub.value = "deep"

    save_path = tmp_path / "complex_node"
    node.to_zarr(save_path)

    # Load and verify
    loaded = open_datanode(save_path)

    assert loaded.num == 123
    assert loaded.text == "hello"
    assert loaded.group.a == 1
    assert loaded.group.b == 2
    assert loaded.group.sub.value == "deep"


def test_cleanup(tmp_path):
    """Clean up test directories"""
    shutil.rmtree(tmp_path, ignore_errors=True)
