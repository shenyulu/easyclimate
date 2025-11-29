easyclimate.core.datanode
=========================

.. py:module:: easyclimate.core.datanode

.. autoapi-nested-parse::

   Hierarchical data structure that dynamically manages attributes and nested nodes



Classes
-------

.. autoapisummary::

   easyclimate.core.datanode.DataNode


Functions
---------

.. autoapisummary::

   easyclimate.core.datanode.open_datanode


Module Contents
---------------

.. py:class:: DataNode(name='root')

   A hierarchical data structure that dynamically manages attributes and nested nodes.

   The :py:class:`DataNode <DataNode>` class provides a flexible way to organize and access data in a tree-like
   structure. It supports automatic creation of nested nodes, path-style access,
   and rich HTML representation in Jupyter environments.

   Parameters
   ----------
   name : :py:class:`str <str>`, optional
       The name of the root node (default: ``"root"``).


   .. py:attribute:: _attributes


   .. py:attribute:: name
      :value: 'root'



   .. py:method:: __getattr__(key)

      Dynamically access or create node attributes.

      Automatically creates nested :py:class:`DataNode <DataNode>` for non-existent attributes,
      while filtering out special IPython/Jupyter attribute requests.

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute name to access.

      Returns
      -------
      Any
          The requested attribute or a new :py:class:`DataNode <DataNode>` if attribute doesn't exist.

      Raises
      ------
      AttributeError
          If the attribute is a special IPython/Jupyter attribute.



   .. py:method:: __setattr__(key, value)

      Set node attributes while protecting internal attributes.

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute name to set.
      value : Any
          The value to assign to the attribute.



   .. py:method:: __getitem__(key)

      Access attributes using path-style notation (e.g., "path/to/attribute").

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute path to access, with '/' separators for nested nodes.

      Returns
      -------
      Any
          The value at the specified path.

      Raises
      ------
      KeyError
          If any part of the path doesn't exist.



   .. py:method:: __contains__(key)


   .. py:method:: __setitem__(key, value)

      Set attributes using path-style notation, creating intermediate nodes as needed.

      Parameters
      ----------
      key : :py:class:`str <str>`
          The attribute path to set, with '/' separators for nested nodes.
      value : Any
          The value to assign at the specified path.



   .. py:method:: _repr_html_()

      Generate HTML representation for Jupyter notebooks.

      Returns
      -------
      :py:class:`str <str>`
          HTML string representing the node and its contents.



   .. py:method:: _format_html()

      Generate complete HTML representation including styles and scripts.

      Returns
      -------
      :py:class:`str <str>`
          Complete HTML document as a string.



   .. py:method:: _format_node_html(node, level=0, parent_id=None)

      Recursively generate HTML for a node and its children.

      Parameters
      ----------
      node : :py:class:`DataNode <DataNode>`
          The node to format.
      level : :py:class:`int <int>`, optional
          Current nesting level (default: 0).
      parent_id : :py:class:`str <str>`, optional
          ID of parent node for DOM construction (default: None).

      Returns
      -------
      :py:class:`str <str>`
          HTML string representing the node.



   .. py:method:: _format_value(value)

      Format values for display, truncating long sequences.

      Parameters
      ----------
      value : Any
          The value to format.

      Returns
      -------
      :py:class:`str <str>`
          Formatted string representation of the value.



   .. py:method:: format_tree(level=0, html=False, is_last_list=None)

      Generate a tree-structured representation of the node.

      Parameters
      ----------
      level : :py:class:`int <int>`, optional
          Current indentation level (default: 0).
      html : :py:class:`bool <bool>`, optional
          Whether to generate HTML output (default: False).
      is_last_list : :py:class:`list <list>` of :py:class:`bool <bool>`, optional
          Track position in hierarchy for proper indentation (default: None).

      Returns
      -------
      :py:class:`str <str>`
          Formatted tree representation.



   .. py:method:: __repr__()


   .. py:method:: _repr_pretty_(p, cycle)

      Support the IPython of pretty printing



   .. py:method:: _repr_mimebundle_(include=None, exclude=None)


   .. py:method:: __dir__()


   .. py:method:: to_zarr(filepath: Union[str, pathlib.Path], zarr_format: Literal[2, 3] = 2)

      Save the DataNode and its contents to a Zarr storage format.

      Parameters
      ----------
      filepath : Union[str, Path]
          Directory path to save the data.
      zarr_format : Literal[2, 3], optional
          Zarr storage format version (default: 2).



   .. py:method:: load(filepath: Union[str, pathlib.Path])
      :classmethod:


      Load a DataNode from a Zarr storage directory.

      Parameters
      ----------
      filepath : Union[str, Path]
          Directory path containing the saved DataNode.

      Returns
      -------
      :py:class:`DataNode <DataNode>`
          The reconstructed DataNode with all its attributes.



.. py:function:: open_datanode(filepath: str) -> DataNode

   Load a DataNode object from native location.

   This function provides a convenient way to load a DataNode that was previously saved
   using the ``DataNode.to_zarr()`` method.

   Parameters
   ----------
   filepath : :py:class:`str <str>`
       The path to the directory containing the saved DataNode data.
       This should be the same path used with DataNode.to_zarr().

   Returns
   -------
   :py:class:`DataNode <DataNode>`
       The loaded DataNode object with all its attributes and nested structure.

   Examples
   --------
   >>> node = open_datanode("path/to/saved_node")
   >>> node.some_attribute  # Access attributes as usual

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_multieof.py


