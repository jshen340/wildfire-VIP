#include <iostream>
#include <fstream>
#include "Field.h"
#include "Mesh.h"
#include "ArrayIO.h"
#include "PolyscopeViewer.h"
#include "PolyscopeViewerColor.h"

void PolyscopeViewer::Initialize(const json& j) {
    // Set up Polyscope viewer options
    const json& j_settings = j.at("settings");
    // TODO: set up settins
    polyscope::view::upDir = polyscope::view::UpDir::YUp;
    polyscope::view::bgColor = { 0.07, 0.07, 0.07, 0.0 };
    polyscope::options::programName = "PolypleX";
    polyscope::options::printPrefix = "[PolypleX] ";
    polyscope::options::autocenterStructures = false;
    polyscope::options::autoscaleStructures = false;
    polyscope::options::groundPlaneEnabled = false;
    polyscope::options::maxFPS = 30;
    polyscope::options::ssaaFactor = 3;

    // Initialize polyscope
    polyscope::init();

    // Read meta data
    Register_Common_Data();
    j_fields = j.at("fields");
    Register_Data();

    // Specify the callback
    polyscope::state::userCallback = std::bind(&PolyscopeViewer::Run, this);
}

void PolyscopeViewer::Register_Common_Data() {
    Update_Last_Frame();
    if (verbose) std::cout << "Last frame: " << last_frame << std::endl;
}

void PolyscopeViewer::Register_Data() {
    for (auto it = j_fields.begin(); it != j_fields.end(); it++) {
        const std::string& field_name = it.key();
        const json& j_field_settings = it.value();
        std::string type = j_field_settings.value<std::string>("type", "");
        bool enabled = j_field_settings.value<bool>("enabled", false);
        if (type == "grid scalar") {
            std::string grid_name = j_field_settings.value<std::string>("grid name", "grid");
            std::string color_map = j_field_settings.value<std::string>("color map", "jet");
            Register_Grid_Cell_Scalar_Field(grid_name, field_name, enabled, color_map);
        }
        else if (type == "grid vector") {
            std::string grid_name = j_field_settings.value<std::string>("grid name", "grid");
            Register_Grid_Cell_Vector_Field(grid_name, field_name, enabled, PolyscopeViewerColor::DeepBlue);
        }
        else if (type == "point cloud") {
            Register_Point_Cloud(field_name, enabled, PolyscopeViewerColor::DeepBlue);
        }
        else {
            std::cout << "Undefined field " << field_name << " with type " << type << std::endl;
        }
    }
}

void PolyscopeViewer::Update_Data() {
    Update_Last_Frame();
    for (auto it = j_fields.begin(); it != j_fields.end(); it++) {
        const std::string& field_name = it.key();
        const json& j_field_settings = it.value();
        std::string type = j_field_settings.value<std::string>("type", "");
        if (type == "grid scalar") {
            const std::string& grid_name = j_field_settings.value<std::string>("grid name", "grid");
            Update_Grid_Cell_Scalar_Field(grid_name, field_name);
        }
        else if (type == "grid vector") {
            const std::string& grid_name = j_field_settings.value<std::string>("grid name", "grid");
            Update_Grid_Cell_Vector_Field(grid_name, field_name);
        }
        else if (type == "point cloud") {
            Update_Point_Cloud(field_name);
        }
    }
}

void PolyscopeViewer::Register_UI()
{
    ImGui::Text("viewer mode: Base");
}

void PolyscopeViewer::Run() {
    if (play) frame++;
    frame = std::max(frame, 0);
    frame = std::min(frame, last_frame);

    ImGui::PushItemWidth(100);

    // Current frame
    ImGui::InputInt("Frame", &frame);

    // Last frame button
    if (ImGui::Button("<", ImVec2(30, 30))) frame = std::max(frame-1, 0);

    // Play/Pause button
    ImGui::SameLine();
    auto pause_name = play ? "Pause" : "Play";
    if (ImGui::Button(pause_name, ImVec2(60, 30))) play = !play;

    // Next frame button
    ImGui::SameLine();
    if (ImGui::Button(">", ImVec2(30, 30))) frame = std::min(frame+1, last_frame);

    // Reset button
    ImGui::SameLine();
    if (ImGui::Button("Reset", ImVec2(60, 30))) frame = 0;

    // Register user data and UI
    if (loaded_frame != frame) {
        Update_Data();
        if (verbose) std::cout << "Read frame " << frame << std::endl;
        loaded_frame = frame;
    }
    Register_UI();

    // Render
    if (ImGui::Checkbox("Render", &render)) Render();

    ImGui::PopItemWidth();
}

void PolyscopeViewer::Render() {
    // Render
    polyscope::screenshot(render_dir);
}

void PolyscopeViewer::Show() {
    // Display
    polyscope::show();
}

polyscope::PointCloud* PolyscopeViewer::Register_Point_Cloud(std::string name, bool enabled, glm::vec3 color) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        Array<float> points_alpha;
        BinaryDataIO::Read_Array(file_name, points_alpha);
        if (verbose) std::cout << "Read file " << name << std::endl;
        Array<Vector3f> points;
        Extract_Points(points_alpha, points);
        polyscope::PointCloud* point_cloud = polyscope::registerPointCloud(name, points);
        point_cloud->setEnabled(enabled);
        point_cloud->setPointRadius(0.001);
        point_cloud->setPointColor(color);
        point_cloud->setPointRenderMode(polyscope::PointRenderMode::Quad);
        point_cloud->setMaterial("flat");
        if (verbose) std::cout << "Register field " << name << std::endl;
        return point_cloud;
    }
    return nullptr;
}

polyscope::PointCloud* PolyscopeViewer::Update_Point_Cloud(std::string name) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        Array<float> points_alpha;
        BinaryDataIO::Read_Array(file_name, points_alpha);
        if (verbose) std::cout << "Read file " << name << std::endl;
        Array<Vector3f> points;
        Extract_Points(points_alpha, points);
        polyscope::PointCloud* point_cloud = polyscope::registerPointCloud(name, points);
        if (verbose) std::cout << "Update field " << name << std::endl;
        return point_cloud;
    }
    return nullptr;
}

PolyscopeGrid* PolyscopeViewer::Register_Grid(std::string name, bool enabled) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        PolyscopeGrid grid;
        File::Read_Binary_From_File(file_name, grid.mac_grid.grid);
        grid.Initialize(name, grid.mac_grid.grid, enabled);
        if (verbose) std::cout << "Read file " << name << std::endl;
        return &grid;
    }
    return nullptr;
}

void PolyscopeViewer::Register_Grid_Cell_Scalar_Field(std::string grid_name, std::string name, bool enabled, std::string color_map) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        polyscope::PointCloud* vis_cell = polyscope::getPointCloud(grid_name + "_cells");
        Field<real, 3> field;
        File::Read_Binary_From_File(file_name, field);
        if (verbose) std::cout << "Read file " << name << std::endl;
        polyscope::PointCloudScalarQuantity* quantity = vis_cell->addScalarQuantity(name, field.array);
        quantity->setEnabled(enabled);
        quantity->setColorMap(color_map);
        if (verbose) std::cout << "Add field " << name << " to " << vis_cell->name << std::endl;
    }
}

void PolyscopeViewer::Update_Grid_Cell_Scalar_Field(std::string grid_name, std::string name) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        polyscope::PointCloud* vis_cell = polyscope::getPointCloud(grid_name + "_cells");
        Field<real, 3> field;
        File::Read_Binary_From_File(file_name, field);
        if (verbose) std::cout << "Read file " << name << std::endl;
        polyscope::PointCloudScalarQuantity* quantity = vis_cell->addScalarQuantity(name, field.array);
        if (verbose) std::cout << "Update field " << name << " in " << vis_cell->name << std::endl;
    }
}

void PolyscopeViewer::Register_Grid_Cell_Vector_Field(std::string grid_name, std::string name, bool enabled, glm::vec3 color) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        polyscope::PointCloud* vis_cell = polyscope::getPointCloud(grid_name + "_cells");
        Field<Vector3, 3> field;
        File::Read_Binary_From_File(file_name, field);
        if (verbose) std::cout << "Read file " << name << std::endl;
        polyscope::PointCloudVectorQuantity* quantity = vis_cell->addVectorQuantity(name, field.array);
        quantity->setEnabled(enabled);
        quantity->setVectorColor(color);
        quantity->setVectorRadius(0.001);
        quantity->setMaterial("flat");
        if (verbose) std::cout << "Add field " << name << " to " << vis_cell->name << std::endl;
    }
}

void PolyscopeViewer::Update_Grid_Cell_Vector_Field(std::string grid_name, std::string name) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        polyscope::PointCloud* vis_cell = polyscope::getPointCloud(grid_name + "_cells");
        Field<Vector3, 3> field;
        File::Read_Binary_From_File(file_name, field);
        if (verbose) std::cout << "Read file " << name << std::endl;
        vis_cell->addVectorQuantity(name, field.array);
        if (verbose) std::cout << "Update field " << name << " in " << vis_cell->name << std::endl;
    }
}

polyscope::CurveNetwork* PolyscopeViewer::Register_Surface_Mesh_2D(std::string name, bool enabled, glm::vec3 color) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        SegmentMesh<3> mesh;
        File::Read_Binary_From_File(file_name, mesh);
        if (verbose) std::cout << "Read file " << name << std::endl;
        polyscope::CurveNetwork* segment_mesh = polyscope::registerCurveNetwork(name, *mesh.vertices, mesh.elements);
        segment_mesh->setColor(color);
        segment_mesh->setMaterial("flat");
        segment_mesh->setRadius(0.001);
        segment_mesh->setEnabled(enabled);
        if (verbose) std::cout << "Register field " << name << std::endl;
        return segment_mesh;
    }
    return nullptr;
}

polyscope::CurveNetwork* PolyscopeViewer::Update_Surface_Mesh_2D(std::string name) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        SegmentMesh<3> mesh;
        File::Read_Binary_From_File(file_name, mesh);
        if (verbose) std::cout << "Read file " << name << std::endl;
        polyscope::CurveNetwork* segment_mesh = polyscope::registerCurveNetwork(name, *mesh.vertices, mesh.elements);
        if (verbose) std::cout << "Update field " << name << std::endl;
        return segment_mesh;
    }
    return nullptr;
}

polyscope::SurfaceMesh* PolyscopeViewer::Register_Surface_Mesh(std::string name, bool enabled, glm::vec3 color) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        TriangleMesh<3> mesh;
        File::Read_Binary_From_File(file_name, mesh);
        if (verbose) std::cout << "Read file " << name << std::endl;
        polyscope::SurfaceMesh* surface_mesh = polyscope::registerSurfaceMesh(name, *mesh.vertices, mesh.elements);
        surface_mesh->setSurfaceColor(color);
        surface_mesh->setEnabled(enabled);
        if (verbose) std::cout << "Register field " << name << std::endl;
        return surface_mesh;
    }
    return nullptr;
}

polyscope::SurfaceMesh* PolyscopeViewer::Update_Surface_Mesh(std::string name) {
    std::string file_name = fmt::format("{}/{}/{}", output_dir, frame, name);
    if (File::File_Exists(file_name)) {
        TriangleMesh<3> mesh;
        File::Read_Binary_From_File(file_name, mesh);
        if (verbose) std::cout << "Read file " << name << std::endl;
        polyscope::SurfaceMesh* surface_mesh = polyscope::registerSurfaceMesh(name, *mesh.vertices, mesh.elements);
        if (verbose) std::cout << "Register field " << name << std::endl;
        return surface_mesh;
    }
    return nullptr;
}

//////////////////////////////////////////////////////////////////////////
////Helper functions

void PolyscopeViewer::Extract_Points(Array<float>& gl_array, Array<Vector3f>& points) {
    const int n = gl_array.size()/4;
    points.resize(n);
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        points[i] = Vector3f(gl_array[4*i], gl_array[4*i+1], gl_array[4*i+2]);
    }
}

void PolyscopeViewer::Update_Last_Frame() {
    std::ifstream in(output_dir + "/0/last_frame.txt");
    Assert((bool)in, "Cannot find last_frame.txt");
    in >> last_frame;
    in.close();
}
