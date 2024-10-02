//#####################################################################
// Main
// SimpleX, Bo Zhu
//#####################################################################
#include <iostream>
#include "Common.h"
//#include "AuxFunc.h"
//#include "Hashtable.h"
//#include "Grid.h"
//#include "MacGrid.h"
//#include "Field.h"
//#include "FaceField.h"
//#include "Interpolation.h"
//#include "Particles.h"
//#include "Mesh.h"
//#include "File.h"
//#include "SparseFunc.h"
//#include "TinyObjLoader.h"
#ifdef USE_CUDA
#include "test_cuda.h"
#endif

//void Test_Misc()
//{
//    Vector3 v(1,2,3);
//    Vector2 v2(1,2);
//    std::cout<<v.transpose()<<", "<<v2.transpose()<<std::endl;
//
//    Vector3d vd(1.,2.,3.);
//    Vector3i vi=vd.cast<int>();
//    Vector3d vt=vi.cast<real>();
//    std::cout<<"vi: "<<vi.transpose()<<", vt: "<<vt.transpose()<<std::endl;
//
//    real array[5]={1,2,3,4,5};
//    for(auto x:array)std::cout<<x<<", ";std::cout<<std::endl;
//}
//
//void Test_Hashtable()
//{
//    Hashtable<int,int> hashtable;
//}
//
//void Test_Grid()
//{
//    const int d=3;
//    Typedef_VectorDii(d);
//    VectorD domain_size=VectorD::Ones();VectorDi cell_counts=VectorDi::Ones()*4;real dx=domain_size[0]/(real)(cell_counts[0]);
//
//    std::cout<<"#test grid"<<std::endl;
//    Grid<d> grid(cell_counts,dx);
//    std::cout<<"cell_counts: "<<grid.cell_counts.transpose()<<", dx: "<<grid.dx<<std::endl;
//    std::cout<<"#test grid operations"<<std::endl;
//    int cell_index=10;
//    VectorDi cell=grid.Cell_Coord(cell_index);
//    std::cout<<"coord: "<<cell.transpose()<<std::endl;
//    VectorD center=grid.Center(cell);
//    std::cout<<"center: "<<center.transpose()<<std::endl;
//    std::cout<<"nb_r: "<<std::endl;
//    for(int i=0;i<Grid<d>::Number_Of_Nb_R();i++)std::cout<<"nb "<<i<<": "<<grid.Nb_R(cell,i).transpose()<<std::endl;
//    std::cout<<"cell_inc_nodes: "<<std::endl;
//    for(int i=0;i<Grid<d>::Number_Of_Cell_Incident_Nodes();i++)std::cout<<"node "<<i<<": "<<grid.Cell_Incident_Node(cell,i).transpose()<<std::endl;
//
//    std::cout<<"#test iter cell"<<std::endl;
//    iterate_cell(iter,grid){
//        std::cout<<"cell "<<iter.Index()<<", "<<iter.Coord().transpose()<<std::endl;}
//    std::cout<<"#test iter range"<<std::endl;
//    iterate_range(iter,VectorDi::Ones()*4,VectorDi::Ones()*6){
//        std::cout<<"cell "<<iter.Index()<<", "<<iter.Coord().transpose()<<std::endl;}
//
//    std::cout<<"#test mac grid"<<std::endl;
//    MacGrid<d> mac_grid(cell_counts,dx);
//    iterate_face(axis,iter,mac_grid){
//        std::cout<<"face "<<axis<<", "<<iter.Coord().transpose()<<std::endl;}
//
//    std::cout<<"#test field"<<std::endl;
//    Field<real,d> field(grid.cell_counts);
//    field.Fill((real)1);
//    field+=(real)10;
//    iterate_cell(iter,grid){
//        std::cout<<"field "<<iter.Coord().transpose()<<": "<<field(iter.Coord())<<std::endl;}
//
//    std::cout<<"#test face field"<<std::endl;
//    FaceField<real,d> face_field(mac_grid.grid.cell_counts);
//    face_field.Fill((real)0,0);
//    face_field.Fill((real)1,1);
//    iterate_face(axis,iter,mac_grid){
//        std::cout<<"face "<<axis<<", "<<iter.Coord().transpose()<<": "<<face_field(axis,iter.Coord())<<std::endl;}
//
//    std::cout<<"#test interpolation"<<std::endl;
//    Interpolation<d> intp(mac_grid);
//    VectorD pos=mac_grid.grid.Center(cell);
//    VectorD val=intp.Interpolate_Face_Vectors(face_field,pos);
//    std::cout<<"intp value: "<<val.transpose()<<std::endl;
//}
//
//void Test_Particles_And_Mesh()
//{
//    const int d=3;
//    Typedef_VectorDii(d);
//    Points<3> points;
//    points.Add_Element();
//    points.X(0)=VectorD::Ones();
//    std::cout<<"size: "<<points.Size()<<std::endl;
//    std::cout<<"x0: "<<points.X(0).transpose()<<std::endl;
//
//    Particles<3> particles;
//    particles.Add_Element();
//    particles.X(0)=VectorD::Ones();
//    particles.V(0)=VectorD::Ones()*(real)2;
//    std::cout<<"size: "<<particles.Size()<<std::endl;
//    std::cout<<"x0: "<<particles.X(0).transpose()<<", v0: "<<particles.V(0).transpose()<<std::endl;
//
//    Array<VectorD> vtx={VectorD::Zero(),VectorD::Unit(0),VectorD::Unit(1)};
//    auto vtx_ptr=std::make_shared<Array<VectorD> >();*vtx_ptr=vtx;
//    TriangleMesh<3> triangle_mesh(vtx_ptr);
//    triangle_mesh.elements.push_back(VectorDi(0,1,2));
//    std::cout<<"tri mesh: "<<triangle_mesh.vertices->size()<<", "<<triangle_mesh.elements.size()<<std::endl;
//}
//
//void Test_IO()
//{
//    const int d=3;
//    Typedef_VectorDii(d);
//    VectorD domain_size=VectorD::Ones();VectorDi cell_counts=VectorDi::Ones()*4;real dx=domain_size[0]/(real)(cell_counts[0]);
//
//    std::cout<<"#test IO"<<std::endl;
//    Grid<d> grid(cell_counts,dx);
//
//    bool has_wrt=File::CheckIOFunc<Grid<d> >::has_binary_read;
//    std::cout<<"check "<<has_wrt<<std::endl;
//
//    std::string file_name="grid";
//    File::Write_Binary_To_File<Grid<d> >(file_name,grid);
//    Grid<d> g2;File::Read_Binary_From_File<Grid<d> >(file_name,g2);
//    std::cout<<"g2: "<<g2.cell_counts.transpose()<<", "<<g2.dx<<", "<<g2.domain_min.transpose()<<std::endl;
//}
//
//void Test_Tiny_Obj()
//{
//    Array<std::shared_ptr<TriangleMesh<3> > > meshes;
//    Obj::Read_From_Obj_File(Path::Data()+"/meshes/bunny.obj",meshes);
//}

void Test_Cuda()
{
#ifdef USE_CUDA
	Test();
#endif
}

void Test_Advanced_Json() {
	json j = { {"init_vec",Vector3(1,9,2)} };
	j["added_vec"] = Vector3(6, 0, 8);
	Vector3 v1 = j["init_vec"];
	Vector3 v2 = j.value<Vector3>("added_vec", Vector3::Zero());
	Vector3 v3 = j.value<Vector3>("naive", Vector3::Ones());
	std::vector<int> std_v = { 1,7 };
	fmt::print("print scalar and std::vector: {} {}\n", 42, std_v);
	fmt::print("print Eigen vector v1 (unsafe): {}\n", v1);
	fmt::print("safely retrieve v2: {}\n", v2);
	fmt::print("non-existing key, use default value v3: {}\n", v3);
	j["scalar_value"] = 42;
	fmt::print("scalar value: {}\n", j.value<int>("scalar_value", 0));
	j["string_value"] = "young";
	fmt::print("string value: {}\n", j.value<std::string>("string_value", ""));
}

int main()
{
    //Test_Misc();
    //Test_Grid();
    //Test_Particles_And_Mesh();
    //Test_IO();
    //Test_Tiny_Obj();
	Test_Cuda();
	Test_Advanced_Json();
	if (true) {
		Assert(false, "assert as python style: {}", "python");
	}
    return 0;
}