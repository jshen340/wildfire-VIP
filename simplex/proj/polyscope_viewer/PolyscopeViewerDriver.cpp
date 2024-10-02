#include "PolyscopeViewerDriver.h"
#include "PolyscopeViewerColor.h"

//////////////////////////////////////////////////////////////////////////
////fluid viewers
void PolyscopeViewerEulerianFluid::Register_Common_Data()
{
    Base::Register_Common_Data();
    Register_Grid("grid", true);
}

void PolyscopeViewerEulerianFluid::Register_Data()
{
    Register_Grid_Cell_Scalar_Field("grid", "pressure", false, "jet");
    Register_Grid_Cell_Vector_Field("grid", "velocity", true, PolyscopeViewerColor::DeepBlue);
    Register_Surface_Mesh_2D("segment_mesh", true, PolyscopeViewerColor::Cyan);
    Register_Surface_Mesh("triangle_mesh", true, PolyscopeViewerColor::Cyan);
    Base::Register_Data();
}

void PolyscopeViewerEulerianFluid::Update_Data()
{
    Update_Grid_Cell_Scalar_Field("grid", "pressure");
    Update_Grid_Cell_Vector_Field("grid", "velocity");
    Update_Surface_Mesh_2D("segment_mesh");
    Update_Surface_Mesh("triangle_mesh");
    Base::Update_Data();
}

void PolyscopeViewerEulerianFluid::Register_UI()
{
    ImGui::Text("viewer mode: Eulerian Fluid");
}
