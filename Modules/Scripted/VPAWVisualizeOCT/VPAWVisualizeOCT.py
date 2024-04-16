import logging
import os
from pathlib import Path

import ctk
import numpy as np
import slicer
import slicer.ScriptedLoadableModule
import slicer.util
import vtk

#
# VPAWVisualizeOCT
#


class VPAWVisualizeOCT(slicer.ScriptedLoadableModule.ScriptedLoadableModule):
    """
    Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        slicer.ScriptedLoadableModule.ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "VPAW Visualize OCT"
        self.parent.categories = ["VPAW"]
        # TODO: add here list of module names that this module requires
        self.parent.dependencies = []
        self.parent.contributors = [
            "Lee Newberg (Kitware, Inc.)",
            "Ebrahim Ebrahim (Kitware, Inc.)",
            "Andinet Enquobahrie (Kitware, Inc.)",
        ]
        # TODO: update with short description of the module and a link to online module
        # documentation
        self.parent.helpText = """
This is the scripted loadable module named VPAW Visualize OCT.  See more information in
<a href="https://github.com/KitwareMedical/vpaw#VPAWVisualizeOCT">module
documentation</a>.
"""
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = """
This file was built from template originally developed by Jean-Christophe Fillion-Robin,
Kitware Inc., Andras Lasso, PerkLab, and Steve Pieper, Isomics, Inc. and was partially
funded by NIH grant 3P41RR013218-12S1.
"""
        # Additional initialization step after application startup is complete
        slicer.app.connect("startupCompleted()", registerSampleData)


#
# Register sample data sets in Sample Data module
#


def registerSampleData():
    """
    Add data sets to Sample Data module.
    """
    # It is always recommended to provide sample data for users to make it easy to try
    # the module, but if no sample data are available then this method (and associated
    # startupCompeted signal connection) can be removed.
    pass


#
# VPAWVisualizeOCTWidget
#


class VPAWVisualizeOCTWidget(
    slicer.ScriptedLoadableModule.ScriptedLoadableModuleWidget,
    slicer.util.VTKObservationMixin,
):
    """
    Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent=None):
        """
        Called when the user opens the module the first time and the widget is
        initialized.
        """
        slicer.ScriptedLoadableModule.ScriptedLoadableModuleWidget.__init__(
            self, parent,
        )
        # needed for parameter node observation:
        slicer.util.VTKObservationMixin.__init__(self)
        self.logic = None
        self._parameterNode = None
        self._updatingGUIFromParameterNode = False

    def setup(self):
        """
        Called when the user opens the module the first time and the widget is
        initialized.
        """
        slicer.ScriptedLoadableModule.ScriptedLoadableModuleWidget.setup(self)

        # Load widget from .ui file (created by Qt Designer).  Additional widgets can be
        # instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/VPAWVisualizeOCT.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets.  Make sure that in Qt designer the top-level
        # qMRMLWidget's "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each
        # MRML widget's "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Configure 3D view
        viewNode = slicer.app.layoutManager().threeDWidget(0).mrmlViewNode()
        viewNode.SetBackgroundColor(0, 0, 0)
        viewNode.SetBackgroundColor2(0, 0, 0)
        viewNode.SetAxisLabelsVisible(False)
        viewNode.SetBoxVisible(False)
        viewNode.SetOrientationMarkerType(
            slicer.vtkMRMLAbstractViewNode.OrientationMarkerTypeAxes,
        )

        # Create logic class.  Logic implements all computations that should be possible
        # to run in batch mode, without a graphical user interface.
        self.logic = VPAWVisualizeOCTLogic()

        # Connections

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(
            slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose,
        )
        self.addObserver(
            slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose,
        )

        # These connections ensure that whenever user changes some settings on the GUI,
        # that is saved in the MRML scene (in the selected parameter node).
        self.ui.DataDirectory.connect(
            "currentPathChanged(const QString&)", self.updateParameterNodeFromGUI,
        )
        self.ui.ImageId.connect(
            "textChanged(const QString&)", self.updateParameterNodeFromGUI,
        )
        self.ui.DataDirectory.connect(
            "validInputChanged(bool)", self.updateParameterNodeFromGUI,
        )

        # Buttons
        self.ui.VPAWModelButton.connect("clicked(bool)", self.onVPAWModelButton)
        self.ui.VPAWVisualizeButton.connect("clicked(bool)", self.onVPAWVisualizeButton)
        self.ui.VPAWModelOCTButton.connect("clicked(bool)", self.onVPAWModelOCTButton)
        self.ui.HomeButton.connect("clicked(bool)", self.onHomeButton)
        self.ui.showButton.connect("clicked(bool)", self.onShowButton)

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

        self.ui.subjectHierarchyTree.setMRMLScene(slicer.mrmlScene)
        shNode = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(
            slicer.mrmlScene,
        )
        self.addObserver(
            shNode, shNode.SubjectHierarchyItemModifiedEvent, self.shItemModifiedEvent,
        )

    def cleanup(self):
        """
        Called when the application closes and the module widget is destroyed.
        """
        self.removeObservers()

    def enter(self):
        """
        Called each time the user opens this module.
        """
        # Make sure parameter node exists and observed
        self.initializeParameterNode()

    def exit(self):
        """
        Called each time the user opens a different module.
        """
        # Do not react to parameter node changes (GUI wlil be updated when the user
        # enters into the module)
        self.removeObserver(
            self._parameterNode,
            vtk.vtkCommand.ModifiedEvent,
            self.updateGUIFromParameterNode,
        )

    def onSceneStartClose(self, caller, event):
        """
        Called just before the scene is closed.
        """
        # Parameter node will be reset, do not use it anymore
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event):
        """
        Called just after the scene is closed.
        """
        # If this module is shown while the scene is closed then recreate a new
        # parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()

    def initializeParameterNode(self):
        """
        Ensure parameter node exists and observed.
        """
        # Parameter node stores all user choices in parameter values, node selections,
        # etc. so that when the scene is saved and reloaded, these settings are
        # restored.
        self.setParameterNode(self.logic.getParameterNode())

    def setParameterNode(self, inputParameterNode):
        """
        Set and observe parameter node.  Observation is needed because when the
        parameter node is changed then the GUI must be updated immediately.
        """
        if inputParameterNode:
            self.logic.setDefaultParameters(inputParameterNode)

        # Unobserve previously selected parameter node and add an observer to the newly
        # selected.  Changes of parameter node are observed so that whenever parameters
        # are changed by a script or any other module those are reflected immediately in
        # the GUI.
        if self._parameterNode is not None and self.hasObserver(
            self._parameterNode,
            vtk.vtkCommand.ModifiedEvent,
            self.updateGUIFromParameterNode,
        ):
            self.removeObserver(
                self._parameterNode,
                vtk.vtkCommand.ModifiedEvent,
                self.updateGUIFromParameterNode,
            )
        self._parameterNode = inputParameterNode
        if self._parameterNode is not None:
            self.addObserver(
                self._parameterNode,
                vtk.vtkCommand.ModifiedEvent,
                self.updateGUIFromParameterNode,
            )

        # Initial GUI update
        self.updateGUIFromParameterNode()

    def updateGUIFromParameterNode(self, caller=None, event=None):
        """
        This method is called whenever parameter node is changed.  The module GUI is
        updated to show the current state of the parameter node.
        """
        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause
        # infinite loop)
        self._updatingGUIFromParameterNode = True

        # Update node selectors and sliders
        self.ui.DataDirectory.currentPath = self._parameterNode.GetParameter(
            "DataDirectory",
        )
        self.ui.ImageId.text = self._parameterNode.GetParameter("ImageId")

        # Update buttons states and tooltips
        if os.path.isdir(self.ui.DataDirectory.currentPath):
            # Enable show button
            self.ui.showButton.toolTip = (
                f"Show files from {self.ui.DataDirectory.currentPath!r}"
                + (
                    f" with id {self.ui.ImageId.text!r}"
                    if self.ui.ImageId.text != ""
                    else ""
                )
            )
            self.ui.showButton.enabled = True
        else:
            # Disable show button
            self.ui.showButton.toolTip = (
                "Show is disabled; first select a valid data directory."
            )
            self.ui.showButton.enabled = False

        # All the GUI updates are done
        self._updatingGUIFromParameterNode = False

    def updateParameterNodeFromGUI(self, caller=None, event=None):
        """
        This method is called when the user makes any change in the GUI.  The changes
        are saved into the parameter node (so that they are restored when the scene is
        saved and loaded).
        """
        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        # Modify all properties in a single batch
        wasModified = self._parameterNode.StartModify()
        self._parameterNode.SetParameter(
            "DataDirectory", self.ui.DataDirectory.currentPath,
        )
        self._parameterNode.SetParameter("ImageId", self.ui.ImageId.text)
        self._parameterNode.EndModify(wasModified)

    @vtk.calldata_type(vtk.VTK_LONG)
    def shItemModifiedEvent(self, caller, eventId, callData):
        """
        Callback for when a subject hierarchy item is modified.
        """
        pass

    def onHomeButton(self):
        """
        Switch to the "Home" module when the user clicks the button.
        """
        slicer.util.selectModule("Home")

    def onVPAWModelButton(self):
        """
        Switch to the "VPAW Model" module when the user clicks the button.
        """
        slicer.util.selectModule("VPAWModel")

    def onVPAWVisualizeButton(self):
        """
        Switch to the "VPAW Visualize" module when the user clicks the button.
        """
        slicer.util.selectModule("VPAWVisualize")

    def onVPAWModelOCTButton(self):
        """
        Switch to the "VPAW Model OCT" module when the user clicks the button.
        """
        slicer.util.selectModule("VPAWModelOCT")

    def onShowButton(self):
        """
        When the user clicks the Show button, find the requested files and load them in
        to 3D Slicer's subject hierarchy.
        """
        with slicer.util.tryWithErrorDisplay(
            "Failed to show image data.", waitCursor=True,
        ):
            list_of_files = self.logic.find_and_sort_files_with_id(
                self.ui.DataDirectory.currentPath, self.ui.ImageId.text,
            )
            if len(list_of_files) == 0:
                raise FileNotFoundError(
                    "No images found"
                    + (
                        f" with image id {self.ui.ImageId.text!r}."
                        if self.ui.ImageId.text != ""
                        else "."
                    ),
                )
            self.logic.clearSubject()
            self.logic.loadNodesToSubjectHierarchy(list_of_files, self.ui.ImageId.text)
            self.logic.arrangeView()


#
# VPAWVisualizeOCTLogic
#


class VPAWVisualizeOCTLogic(slicer.ScriptedLoadableModule.ScriptedLoadableModuleLogic):
    """
    This class should implement all the actual computation done by your module.  The
    interface should be such that other python code can import this class and make use
    of the functionality without requiring an instance of the Widget.  Uses
    ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):
        """
        Called when the logic class is instantiated.  Can be used for initializing
        member variables.
        """
        slicer.ScriptedLoadableModule.ScriptedLoadableModuleLogic.__init__(self)
        self.clearSubject()

    def setDefaultParameters(self, parameterNode):
        """
        Initialize parameter node with default settings.
        """
        if not parameterNode.GetParameter("Threshold"):
            parameterNode.SetParameter("Threshold", "100.0")
        if not parameterNode.GetParameter("Invert"):
            parameterNode.SetParameter("Invert", "false")

    def find_files_with_id(self, path, id, include_subjectless=False):
        """
        Find all file names within `path` recursively that contain `id`

        Parameters
        ----------
        path : str
            Initially, the top-level directory to be scanned for files.  For recursive
            calls, it will be a directory or file within the top-level directory's
            hierarchy.
        id: str
            A value such as "1000_" will find all proper files that have basenames that
            contain that string.  If id=="" then all files regardless of name will be
            reported.
        include_subjectless: bool
            If set to True then files not associated with any image will also be
            included.

        Returns
        -------
        When `path` is a file, returns the one-tuple list `[(path, mtime)]` if the
            `path` basename contains `id`; otherwise returns an empty list.
            "mtime" is the modification time returned by os.path.getmtime(path).
        When `path` is a directory, returns the concatenation of the lists generated by
            a recursive call to find_files_with_id for each entry in the directory.
        """
        if os.path.isdir(path):
            # This `path` is a directory.  Recurse to all files and directories within
            # `path` and flatten the responses into a single list.
            response = [
                record
                for sub in os.listdir(path)
                for record in self.find_files_with_id(
                    os.path.join(path, sub), id, include_subjectless,
                )
            ]
        else:
            # This path is not a directory.  Return a list containting the path if the
            # path basename contains `id`, otherwise return an empty list.
            response = [
                (p, os.path.getmtime(p))
                for p in (path,)
                if id in os.path.basename(p)
                or (
                    include_subjectless
                    and (
                        "mean_landmarks" in p
                        or "FilteredControlBlindingLogUniqueScanFiltered" in p
                        or "weighted_perc" in p
                    )
                )
            ]
        return response

    def find_and_sort_files_with_id(self, dataDirectory, imageId):
        """
        Find all file names within `path` recursively that contain `id`, and sort them
        by their modification times

        Parameters
        ----------
        dataDirectory : str
            The top-level directory to be scanned for files.
        imageId: str
            A value such as "1000_" will find all proper files that have basenames that
            contain that string.  If id=="" then all files regardless of name will be
            reported.

        Returns
        -------
        When `path` is a directory, returns a list of the requested files, sorted
        chronologically by their modification times in
        """
        if not (isinstance(dataDirectory, str) and os.path.isdir(dataDirectory)):
            raise ValueError(f"Data directory (value={dataDirectory!r}) is not valid")
        if not (imageId is None or isinstance(imageId, str)):
            raise ValueError(f"Image id (value={imageId!r}) is not valid")
        if imageId is None:
            imageId = ""

        import time

        startTime = time.time()
        logging.info("Processing started")

        list_of_records = self.find_files_with_id(
            dataDirectory, imageId, include_subjectless=False,
        )
        # Sort by modification time
        list_of_records.sort(key=lambda record: record[1])
        # Remove modification times
        list_of_files = [record[0] for record in list_of_records]

        stopTime = time.time()
        logging.info(f"Processing completed in {stopTime-startTime:.2f} seconds")

        return list_of_files

    def loadCenterlineFromP3FileContents(self, contents):
        """
        Load a centerline using the data object written into a P3 file named
        "####_CENTERLINE.p3"

        Parameters
        ----------
        contents : a pair of arrays (centerline_points, centerline_normals).  Currently
            we only use centerline_points, piecing them together into a curve node.

        Returns
        -------
        A vtkMRMLMarkupsCurveNode
        """
        centerline_points, centerline_normals = contents

        # The axis ordering is not IJK to begin with, hence this permuation
        centerline_points = centerline_points[:, [2, 1, 0]]

        centerline_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsCurveNode")
        centerline_node.SetName("Centerline")
        slicer.util.updateMarkupsControlPointsFromArray(
            centerline_node, centerline_points,
        )
        centerline_node.GetDisplayNode().SetGlyphTypeFromString("Vertex2D")
        centerline_node.SetCurveTypeToLinear()
        # don't allow mouse interaction to move control points
        centerline_node.LockedOn()
        # hide the text label because it distracts from landmarks
        centerline_node.GetDisplayNode().SetPropertiesLabelVisibility(False)
        return centerline_node

    def loadOneNode(self, filename, basename_repr, props):
        """
        Create a 3D Slicer node object for the data in a file

        Parameters
        ----------
        filename : str
            The file from which to read the data that will define the created node
        basename_repr: str
            A string representing the file, which is used for warning/error output
        props : dict
            A dictionary of properties that is passed to most slicer.util.load*
            functions

        Returns
        -------
        A 3D Slicer node object representing the data
        """
        # Determine the node type from the filename extension, using its
        # immediate-ancestor directory's name if necessary.  Note: check for ".seg.nrrd"
        # before checking for ".nrrd".
        if filename.endswith(".seg.nrrd"):
            node = slicer.util.loadSegmentation(filename, properties=props)
            node.CreateClosedSurfaceRepresentation()
        elif filename.endswith(".nrrd"):
            directory = os.path.basename(os.path.dirname(filename))
            if directory == "images":
                node = slicer.util.loadVolume(filename, properties=props)
                self.show_nodes.append(node)
            elif directory == "segmentations_computed":
                node = slicer.util.loadSegmentation(filename, properties=props)
                node.CreateClosedSurfaceRepresentation()
            else:
                # Guess
                node = slicer.util.loadVolume(filename, properties=props)
        elif filename.endswith(".fcsv"):
            node = slicer.util.loadMarkups(filename)
            assert node.IsTypeOf("vtkMRMLMarkupsNode")
            node.LockedOn()  # don't allow mouse interaction to move control points
        elif filename.endswith(".mha") or filename.endswith(".png"):
            node = slicer.util.loadVolume(filename, properties=props)
        elif filename.endswith(".xls"):
            print(f"File type for {basename_repr} is not currently supported")
            node = None
        else:
            print(f"File type for {basename_repr} is not recognized")
            node = None
        return node

    def clearSubject(self):
        """
        Set VPAWVisualizeOCTLogic to initial state before any subject was loaded, and
        clear the subject hierarchy.
        """
        self.subject_id = None
        # subject hierarchy item id for the currently loaded subject
        self.subject_item_id = None
        self.input_image_node = None
        self.input_ijk_to_ras = None
        self.centerline_node = None
        self.segmentation_node = None
        self.laplace_sol_node = None
        self.laplace_sol_masked_node = None
        self.clearSubjectHierarchy()

    def clearSubjectHierarchy(self):
        """
        Remove all nodes from the 3D Slicer subject hierarchy
        """
        slicer.mrmlScene.GetSubjectHierarchyNode().RemoveAllItems(True)
        self.show_nodes = list()

    def subjectIsCurrentlyLoaded(self) -> bool:
        """
        Whether a subject has been loaded.
        """
        return self.subject_id is not None

    def put_node_under_subject(self, node):
        """
        Move the given node in the subject hierarchy such that it becomes a child of the
        subject item for the currently loaded subject.
        """
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        node_item = shNode.GetItemByDataNode(node)
        shNode.SetItemParent(node_item, self.subject_item_id)

    def loadOneNodeToSubjectHierarchy(self, shNode, subject_item, filename):
        """
        Load data from a single file into a node and put the node in the 3D Slicer
        subject hierarchy

        Parameters
        ----------
        shNode: subject hierarchy node
            The 3D subject hierarchy node
        subject_item: int
            Parent for the node we are creating
        filename: str
            The data source for the file
        """
        # The node types supported by 3D Slicer generally can be found with fgrep
        # 'loadNodeFromFile(filename' from
        # https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/util.py.
        # Currently they are AnnotationFile, ColorTableFile, FiberBundleFile,
        # MarkupsFile, ModelFile, ScalarOverlayFile, SegmentationFile,
        # ShaderPropertyFile, TableFile, TextFile, TransformFile, and VolumeFile.

        basename = os.path.basename(filename)
        basename_repr = repr(basename)
        props = {"name": basename, "singleFile": True, "show": False}

        node = self.loadOneNode(filename, basename_repr, props)
        if node is None:
            return

        dirname = Path(filename).parent.stem
        if dirname == "sols":
            self.laplace_sol_node = node
        elif dirname == "images":
            self.input_image_node = node
        elif dirname == "centerline":
            self.centerline_node = node
        elif dirname == "segmentations_computed":
            self.segmentation_node = node

        self.put_node_under_subject(node)

    def loadNodesToSubjectHierarchy(self, list_of_files, subject_name):
        """
        Load data from files into nodes and put the nodes in the 3D Slicer subject
        hierarchy

        Parameters
        ----------
        list_of_files : List[str]
            Files to be loaded
        subject_name : str
            Name for folder in subject hierarchy to contain the nodes
        """
        self.subject_id = subject_name

        # The subject hierarchy node can contain subject (image), study (optionally),
        # and node items.  slicer.mrmlScene knows how to find the subject hierarchy
        # node.
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        # A subject item is created with the subject hierarchy node as its parent.
        self.subject_item_id = shNode.CreateSubjectItem(
            shNode.GetSceneItemID(), subject_name,
        )

        # slicer knows how to find the subject hierarchy tree view.
        shTV = slicer.qMRMLSubjectHierarchyTreeView()
        # Tell the subject hierarchy tree view about its enclosing scene.
        shTV.setMRMLScene(slicer.mrmlScene)
        # Tell the subject hierarchy tree view that its root item is the subject item.
        shTV.setRootItem(self.subject_item_id)

        for filename in list_of_files:
            self.loadOneNodeToSubjectHierarchy(shNode, self.subject_item_id, filename)

        # further processing that can occur now that all nodes are loaded
        self.create_input_ijk2ras_as_node()
        self.fix_image_origins_and_spacings()
        self.restrict_laplace_sol_to_segmentation()

        # Recursively set visibility and expanded properties of each item
        def recurseVisibility(item, visibility, expanded):
            # Useful functions for traversing items
            # shNode.GetSceneItemID()
            # shNode.GetNumberOfItems()
            # shNode.GetNumberOfItemChildren(parentItem)
            # shNode.GetItemByPositionUnderParent(parentItem, childIndex)
            # shNode.SetItemExpanded(shNode.GetSceneItemID(), True)
            shNode.SetItemDisplayVisibility(item, visibility)
            shNode.SetItemExpanded(item, expanded)
            for child_index in range(shNode.GetNumberOfItemChildren(item)):
                recurseVisibility(
                    shNode.GetItemByPositionUnderParent(item, child_index),
                    visibility,
                    expanded,
                )

        recurseVisibility(self.subject_item_id, True, True)

        # Resize columns of the SubjectHierarchyTreeView
        shTV.header().resizeSections(shTV.header().ResizeToContents)
        # Force re-displaying of the SubjectHierarchyTreeView
        slicer.mrmlScene.StartState(slicer.vtkMRMLScene.ImportState)
        slicer.mrmlScene.EndState(slicer.vtkMRMLScene.ImportState)

    def create_input_ijk2ras_as_node(self):
        """
        Get the IJK to RAS matrix for the input image as a transform node.
        """
        if self.input_image_node is None:
            raise RuntimeError("Could not find input image node.")
        ijkToRas = vtk.vtkMatrix4x4()
        self.input_image_node.GetIJKToRASMatrix(ijkToRas)
        ijkToRas_node = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLLinearTransformNode")
        ijkToRas_node.SetName(f"{self.input_image_node.GetName()}_IJK_to_RAS")
        ijkToRas_node.SetMatrixTransformToParent(ijkToRas)
        self.put_node_under_subject(ijkToRas_node)
        self.input_ijk_to_ras = ijkToRas_node

    def fix_image_origins_and_spacings(self):
        """
        Some nodes rely on others for origin and spacing info, because it wasn't
        properly saved in the files from which we generate those nodes.  This functions
        goes through and transfers origin and spacing info whereever it is needed.
        """
        if self.laplace_sol_node is None:
            raise RuntimeError("Could not find laplace solution node.")
        if self.input_image_node is None:
            raise RuntimeError("Could not find input image node.")
        if self.input_ijk_to_ras is None:
            raise RuntimeError("IJK to RAS transform node has not been created.")
        if self.centerline_node is None:
            raise RuntimeError("Could not find centerline node.")

        self.laplace_sol_node.SetOrigin(self.input_image_node.GetOrigin())
        self.laplace_sol_node.SetSpacing(self.input_image_node.GetSpacing())

        self.centerline_node.SetAndObserveTransformNodeID(self.input_ijk_to_ras.GetID())

    def restrict_laplace_sol_to_segmentation(self):
        """
        If the laplace solution and the segmentation node both exist, mask the laplace
        solution volume by the segmentation node.  If either of them doesn't exists,
        raise an exception.
        """
        if self.segmentation_node is None:
            raise RuntimeError("Could not find segmentation node.")
        if self.laplace_sol_node is None:
            raise RuntimeError("Could not find laplace solution node.")

        sol_array = slicer.util.arrayFromVolume(self.laplace_sol_node)

        seg_ids = self.segmentation_node.GetSegmentation().GetSegmentIDs()
        if len(seg_ids) != 1:
            raise RuntimeError(
                f"Expected node {self.segmentation_node.GetName()} to have"
                " exactly one segmentation.",
            )
        seg_array = slicer.util.arrayFromSegmentBinaryLabelmap(
            self.segmentation_node, seg_ids[0], self.laplace_sol_node,
        )

        sol_masked_array = np.copy(sol_array)
        sol_masked_array[seg_array == 0] = np.nan

        sol_masked_node = slicer.modules.volumes.logic().CloneVolume(
            self.laplace_sol_node,
            self.laplace_sol_node.GetName() + "_restrictedToSegmentation",
        )
        slicer.util.updateVolumeFromArray(sol_masked_node, sol_masked_array)
        self.put_node_under_subject(sol_masked_node)
        self.laplace_sol_masked_node = sol_masked_node

    def arrangeView(self):
        """
        Make the 3D Slicer viewing panels default to something reasonable
        """
        # Make sure at least one input image (if any) is being viewed
        if self.show_nodes:
            slicer.util.setSliceViewerLayers(foreground=self.show_nodes[0], fit=True)
        self.show_nodes = list()

        # Center the 3D view
        layoutManager = slicer.app.layoutManager()
        threeDWidget = layoutManager.threeDWidget(0)
        threeDView = threeDWidget.threeDView()
        threeDView.lookFromAxis(ctk.ctkAxesWidget.Left)
        threeDView.resetFocalPoint()
        threeDView.resetCamera(False, False, True)
        threeDView.lookFromAxis(ctk.ctkAxesWidget.Left)


#
# VPAWVisualizeOCTTest
#


class VPAWVisualizeOCTTest(slicer.ScriptedLoadableModule.ScriptedLoadableModuleTest):
    """
    This is the test case for your scripted module.  Uses ScriptedLoadableModuleTest
    base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def setUp(self):
        """
        Do whatever is needed to reset the state - typically a scene clear will be
        enough.
        """
        slicer.mrmlScene.Clear()

    def runTest(self):
        """
        Run as few or as many tests as needed here.
        """
        self.setUp()
        self.test_VPAWVisualizeOCT1()

    def test_VPAWVisualizeOCT1(self):
        """
        Ideally we should have several levels of tests.  At the lowest level tests
        should exercise the functionality of the logic with different inputs (both valid
        and invalid).  At higher levels our tests should emulate the way the user would
        interact with our code and confirm that it still works the way we intended.

        One of the most important features of the tests is that it should alert other
        developers when their changes will have an impact on the behavior of our module.
        For example, if a developer removes a feature that we depend on, our test should
        break so they know that the feature is needed.
        """
        self.delayDisplay("Starting the test")

        # Get/create input data
        self.delayDisplay("Test skipped")
