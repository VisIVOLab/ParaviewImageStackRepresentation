/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile$

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkVLVAImageStackRepresentation.h"

#include "vtkExtractVOI.h"
#include "vtkImageData.h"
#include "vtkImageMapToColors.h"
#include "vtkImageProperty.h"
#include "vtkImageSliceCollection.h"
#include "vtkImageSliceMapper.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "vtkMPIImageReader.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkProperty.h"
#include "vtkPVLODActor.h"
#include "vtkPVRenderView.h"
#include "vtkRenderer.h"
#include "vtkScalarsToColors.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredExtent.h"
#include "vtkPVTransform.h"

#include <algorithm>


vtkStandardNewMacro(vtkVLVAImageStackRepresentation)
//----------------------------------------------------------------------------
vtkVLVAImageStackRepresentation::vtkVLVAImageStackRepresentation()
{
//    this->SetRepresentation(VTK_SURFACE);
    std::cout << "Initialising image stack representation." << std::endl;
    this->Slice = 0;
    this->SliceMode = XY_PLANE;
    this->SliceMapper = vtkVLVAImageStackMapper::New();
    this->Actor = vtkPVLODActor::New();
    this->Actor->SetMapper(this->SliceMapper);
    this->Actor->GetProperty()->LightingOff();

    vtkSmartPointer<vtkImageSliceMapper> imageSliceMapperBase = vtkSmartPointer<vtkImageSliceMapper>::New();

    std::cout << "Creating image stack base image" << std::endl;
    this->ImageStack = vtkSmartPointer<vtkImageStack>::New();
    auto imageSliceBase = vtkSmartPointer<vtkImageSlice>::New();
    imageSliceBase->SetMapper(imageSliceMapperBase);
    imageSliceBase->GetProperty()->SetInterpolationTypeToNearest();
    imageSliceBase->GetProperty()->SetLayerNumber(0);
    this->ImageStack->AddImage(imageSliceBase);

//    AddLayerImage("/home/pietersielie-aw-ubuntu/Documents/INAF/rosat_pspc_rdf2_3_im3.fits");
}

//----------------------------------------------------------------------------
vtkVLVAImageStackRepresentation::~vtkVLVAImageStackRepresentation()
{
    this->SliceMapper->SetInputData(nullptr);
    this->SliceMapper->Delete();
    this->Actor->Delete();
    this->ImageStack->Delete();
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetSliceMode(int mode)
{
    if (this->SliceMode != mode)
    {
        this->SliceMode = mode;
        this->MarkModified();
    }
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::AddLayerImage(const std::string file)
{
    std::cerr << "Adding file " << file << " to stack." << std::endl;
    if (file.empty()){
        std::cerr << "Adding new image aborted due to empty file name." << std::endl;
        return;
    }
    vtkSmartPointer<vtkMPIImageReader> reader = vtkSmartPointer<vtkMPIImageReader>::New();
    reader->SetFileName(file.c_str());
    std::cerr << "Reader created for file to add to stack." << std::endl;

    vtkSmartPointer<vtkImageMapToColors> colors = vtkSmartPointer<vtkImageMapToColors>::New();
    colors->SetInputData(reader->GetOutput());

    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfTableValues(128);
    for (int i = 0; i < 128; i++)
        lut->SetTableValue(i, i / 128.0, i / 128.0, i / 128.0, 1);
    colors->SetLookupTable(lut);
    colors->Update();

    vtkSmartPointer<vtkImageSliceMapper> imageSliceMapperLayer =
        vtkSmartPointer<vtkImageSliceMapper>::New();
    imageSliceMapperLayer->SetInputData(colors->GetOutput());

    vtkSmartPointer<vtkImageSlice> newImage = vtkSmartPointer<vtkImageSlice>::New();
    newImage->SetMapper(imageSliceMapperLayer);
    newImage->GetProperty()->SetOpacity(0.8);
    newImage->GetProperty()->SetInterpolationTypeToNearest();

    double bounds[6];
    double angle = 0;
    vtkSmartPointer<vtkPVTransform> transform = vtkSmartPointer<vtkPVTransform>::New();
    reader->GetOutput()->GetBounds(bounds);

    // Rotate about the origin point (world coordinates)
    transform->Translate(bounds[0], bounds[2], bounds[4]);
    transform->RotateWXYZ(angle, 0, 0, 1);
    transform->Translate(-bounds[0], -bounds[2], -bounds[4]);
    newImage->SetUserTransform(transform);

    int newActiveLayerNumber = this->ImageStack->GetImages()->GetNumberOfItems();
    newImage->GetProperty()->SetLayerNumber(newActiveLayerNumber);
    this->ImageStack->AddImage(newImage);
    this->SetStackActiveLayer(newActiveLayerNumber);
    std::cerr << "Image " << file << " added to stack." << std::endl;

    this->MarkModified();
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::RemoveLayerImage(const int index)
{
    std::cout << "Removing image that is at index " << index << " from image stack." << std::endl;
    if (index > 0 && index < ImageStack->GetImages()->GetNumberOfItems()){
        auto img = vtkImageSlice::SafeDownCast(ImageStack->GetImages()->GetItemAsObject(index));
        this->ImageStack->RemoveImage(img);
        for (int i = index; i < ImageStack->GetImages()->GetNumberOfItems(); ++i){
            vtkImageSlice::SafeDownCast(this->ImageStack->GetImages()->GetItemAsObject(i))->GetProperty()->SetLayerNumber(i - 1);
        }
        std::cerr << "Image that was at index " << index << " removed from image stack." << std::endl;
    }
    else{
        std::cerr << "Index " << index << " for image to be removed is out of bounds [1, " << ImageStack->GetImages()->GetNumberOfItems() - 1 << "]." << std::endl;
    }

    this->MarkModified();
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetStackLookupTable(vtkScalarsToColors *val)
{
    std::cerr << "Stack lookup table set to " << val << "." << std::endl;
    this->ImageStack->GetActiveImage()->GetProperty()->SetLookupTable(val);

    this->MarkModified();
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetStackOpacity(double val)
{
    std::cerr << "Stack active image opacity set to " << val << "." << std::endl;
    this->ImageStack->GetActiveImage()->GetProperty()->SetOpacity(val);

    this->MarkModified();
}

//----------------------------------------------------------------------------
int vtkVLVAImageStackRepresentation::GetStackActiveLayer() const
{
    std::cerr << "Stack active image index is " << this->ImageStack->GetActiveLayer() << "." << std::endl;
    return this->ImageStack->GetActiveLayer();
}

//----------------------------------------------------------------------------
int vtkVLVAImageStackRepresentation::GetStackLayerCount() const
{
    std::cerr << "Returning count of stack layers which is " << this->ImageStack->GetImages()->GetNumberOfItems() << "." << std::endl;
    return this->ImageStack->GetImages()->GetNumberOfItems();
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetStackActiveLayer(int val)
{
    std::cerr << "Stack active image index set to " << val << "." << std::endl;
    if (val > 0 && val < ImageStack->GetImages()->GetNumberOfItems()){
        this->ImageStack->SetActiveLayer(val);

        this->MarkModified();
    }
    else{
        std::cerr << "Index " << val << " for image to be set as the active layer is out of bounds [1, " << ImageStack->GetImages()->GetNumberOfItems() - 1 << "]." << std::endl;
    }
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetStackLayerVisible(int index, int visible)
{
    if (index > 0 && index < ImageStack->GetImages()->GetNumberOfItems()){
        auto img = vtkImageSlice::SafeDownCast(ImageStack->GetImages()->GetItemAsObject(index));
        img->SetVisibility(visible);

        this->MarkModified();
    }
    else{
        std::cerr << "Index " << index << " for image to have visibility changed is out of bounds [1, " << ImageStack->GetImages()->GetNumberOfItems() - 1 << "]." << std::endl;
    }
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetLayerIndex(int index, int newIndex)
{
    if (index <= 0 || index >= ImageStack->GetImages()->GetNumberOfItems() || newIndex <= 0 || newIndex >= ImageStack->GetImages()->GetNumberOfItems()){
        std::cerr << "Index " << index << " for image to have visibility changed to new index " << newIndex << " is out of bounds [1, " << ImageStack->GetImages()->GetNumberOfItems() - 1 << "]." << std::endl;
    }
    else{
        if (index == newIndex)
            return;
        if (index < newIndex){
            for (int i = index + 1; i <= newIndex; ++i){
                vtkImageSlice::SafeDownCast(this->ImageStack->GetImages()->GetItemAsObject(i))->GetProperty()->SetLayerNumber(i - 1);
            }
        }
        else{
            for (int i = newIndex; i < index; ++i){
                vtkImageSlice::SafeDownCast(this->ImageStack->GetImages()->GetItemAsObject(i))->GetProperty()->SetLayerNumber(i + 1);
            }
        }
        vtkImageSlice::SafeDownCast(this->ImageStack->GetImages()->GetItemAsObject(index))->GetProperty()->SetLayerNumber(newIndex);

        this->MarkModified();
    }
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetSlice(unsigned int val)
{
    if (this->Slice != val)
    {
        this->Slice = val;
        this->MarkModified();
    }
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetInputArrayToProcess(
        int idx, int port, int connection, int fieldAssociation, const char* name)
{
    this->Superclass::SetInputArrayToProcess(idx, port, connection, fieldAssociation, name);
    this->SliceMapper->SelectColorArray(name);
    switch (fieldAssociation)
    {
    case vtkDataObject::FIELD_ASSOCIATION_CELLS:
        this->SliceMapper->SetScalarMode(VTK_SCALAR_MODE_USE_CELL_FIELD_DATA);
        break;

    case vtkDataObject::FIELD_ASSOCIATION_NONE:
        this->SliceMapper->SetScalarMode(VTK_SCALAR_MODE_USE_FIELD_DATA);
        // Color entire block by zeroth tuple in the field data
        this->SliceMapper->SetFieldDataTupleId(0);
        break;

    case vtkDataObject::FIELD_ASSOCIATION_POINTS:
    default:
        this->SliceMapper->SetScalarMode(VTK_SCALAR_MODE_USE_POINT_FIELD_DATA);
        break;
    }
}

//----------------------------------------------------------------------------
int vtkVLVAImageStackRepresentation::FillInputPortInformation(int, vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
}

//----------------------------------------------------------------------------
int vtkVLVAImageStackRepresentation::ProcessViewRequest(
        vtkInformationRequestKey* request_type, vtkInformation* inInfo, vtkInformation* outInfo)
{
    if (!this->Superclass::ProcessViewRequest(request_type, inInfo, outInfo))
    {
        return 0;
    }

    if (request_type == vtkPVView::REQUEST_UPDATE())
    {
        // provide the "geometry" to the view so the view can delivery it to the
        // rendering nodes as and when needed.
        vtkPVView::SetPiece(inInfo, this, this->SliceData);
        vtkPVRenderView::SetGeometryBounds(inInfo, this, this->SliceData->GetBounds());

        // BUG #14253: support translucent rendering.
        if (this->Actor->HasTranslucentPolygonalGeometry())
        {
            outInfo->Set(vtkPVRenderView::NEED_ORDERED_COMPOSITING(), 1);
            // Pass on the partitioning information to the view. This logic is similar
            // to what we do in vtkImageVolumeRepresentation.
            vtkPVRenderView::SetOrderedCompositingConfiguration(
                        inInfo, this, vtkPVRenderView::USE_BOUNDS_FOR_REDISTRIBUTION, this->WholeBounds);
        }
    }
    else if (request_type == vtkPVView::REQUEST_RENDER())
    {
        // since there's no direct connection between the mapper and the collector,
        // we don't put an update-suppressor in the pipeline.
        vtkAlgorithmOutput* producerPort = vtkPVRenderView::GetPieceProducer(inInfo, this);
        if (producerPort)
        {
            this->SliceMapper->SetInputConnection(producerPort);
            this->Actor->GetProperty()->ShadingOff();
        }
    }
    return 1;
}

//----------------------------------------------------------------------------
int vtkVLVAImageStackRepresentation::RequestData(
        vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    vtkMath::UninitializeBounds(this->WholeBounds);
    if (inputVector[0]->GetNumberOfInformationObjects() == 1)
    {
        this->UpdateSliceData(inputVector);
    }
    else
    {
        this->SliceData->Initialize();
    }
    return this->Superclass::RequestData(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::UpdateSliceData(vtkInformationVector** inputVector)
{
    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkImageData* input = vtkImageData::GetData(inputVector[0], 0);

    input->GetBounds(this->WholeBounds);

    int inWholeExtent[6], outExt[6];
    memset(outExt, 0, sizeof(int) * 6);

    inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), inWholeExtent);
    int dataDescription = vtkStructuredData::SetExtent(inWholeExtent, outExt);

    if (vtkStructuredData::GetDataDimension(dataDescription) != 3)
    {
        this->SliceData->ShallowCopy(input);
        return;
    }

    int dims[3];
    dims[0] = inWholeExtent[1] - inWholeExtent[0] + 1;
    dims[1] = inWholeExtent[3] - inWholeExtent[2] + 1;
    dims[2] = inWholeExtent[5] - inWholeExtent[4] + 1;

    unsigned int slice = this->Slice;
    switch (this->SliceMode)
    {
    case YZ_PLANE:
        slice = (static_cast<int>(slice) >= dims[0]) ? dims[0] - 1 : slice;
        outExt[0] = outExt[1] = outExt[0] + slice;
        break;

    case XZ_PLANE:
        slice = (static_cast<int>(slice) >= dims[1]) ? dims[1] - 1 : slice;
        outExt[2] = outExt[3] = outExt[2] + slice;
        break;

    case XY_PLANE:
    default:
        slice = (static_cast<int>(slice) >= dims[2]) ? dims[2] - 1 : slice;
        outExt[4] = outExt[5] = outExt[4] + slice;
        break;
    }

    // Now, clamp the extent for the slice to the extent available on this rank.
    vtkNew<vtkStructuredExtent> helper;
    helper->Clamp(outExt, input->GetExtent());
    if (outExt[0] <= outExt[1] && outExt[2] <= outExt[3] && outExt[4] <= outExt[5])
    {
        vtkExtractVOI* voi = vtkExtractVOI::New();
        voi->SetVOI(outExt);
        voi->SetInputData(input);
        voi->Update();

        this->SliceData->ShallowCopy(voi->GetOutput());
        voi->Delete();
    }
    else
    {
        this->SliceData->Initialize();
    }

    // vtkExtractVOI is passing the origin for the individual slice.
    // We want to use the input origin/spacing to compute the bounds.
    this->SliceData->SetOrigin(input->GetOrigin());
}

//----------------------------------------------------------------------------
bool vtkVLVAImageStackRepresentation::AddToView(vtkView* view)
{
    vtkPVRenderView* rview = vtkPVRenderView::SafeDownCast(view);
    if (rview)
    {
//        rview->GetRenderer()->AddActor(this->Actor);
        rview->GetRenderer()->AddViewProp(this->ImageStack);
        return this->Superclass::AddToView(rview);
    }
    return false;
}

//----------------------------------------------------------------------------
bool vtkVLVAImageStackRepresentation::RemoveFromView(vtkView* view)
{
    vtkPVRenderView* rview = vtkPVRenderView::SafeDownCast(view);
    if (rview)
    {
//        rview->GetRenderer()->RemoveActor(this->Actor);
        rview->GetRenderer()->RemoveViewProp(this->ImageStack);
        return this->Superclass::RemoveFromView(rview);
    }
    return false;
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

//****************************************************************************
// Calls forwarded to internal objects.

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetOrientation(double x, double y, double z)
{
    this->Actor->SetOrientation(x, y, z);
    this->ImageStack->SetOrientation(x, y, z);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetOrigin(double x, double y, double z)
{
    this->Actor->SetOrigin(x, y, z);
    this->ImageStack->SetOrigin(x, y, z);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetPickable(int val)
{
    this->Actor->SetPickable(val);
    this->ImageStack->SetPickable(val);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetPosition(double x, double y, double z)
{
    this->Actor->SetPosition(x, y, z);
    this->ImageStack->SetPosition(x, y, z);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetScale(double x, double y, double z)
{
    this->Actor->SetScale(x, y, z);
    this->ImageStack->SetScale(x, y, z);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetVisibility(bool val)
{
    this->Actor->SetVisibility(val ? 1 : 0);
    this->ImageStack->SetVisibility(val ? 1 : 0);
    this->Superclass::SetVisibility(val);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetOpacity(double val)
{
    this->Actor->GetProperty()->SetOpacity(val);
    this->ImageStack->GetProperty()->SetOpacity(val);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetLookupTable(vtkScalarsToColors* val)
{
    this->SliceMapper->SetLookupTable(val);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetMapScalars(int val)
{
    this->SliceMapper->SetColorMode(val ? VTK_COLOR_MODE_MAP_SCALARS : VTK_COLOR_MODE_DIRECT_SCALARS);
}

//----------------------------------------------------------------------------
void vtkVLVAImageStackRepresentation::SetUseXYPlane(int val)
{
    this->SliceMapper->SetUseXYPlane(val);
}
