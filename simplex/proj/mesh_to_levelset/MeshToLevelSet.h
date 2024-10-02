//////////////////////////////////////////////////////////////////////////
// Copyright (c) (2020-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __MeshToLevelSet_h__
#define __MeshToLevelSet_h__
#include "LevelSet.h"
#include "Mesh.h"
#include "MeshFunc.h"
#include "Common.h"
#include "SimplicialPrimitives.h"
#include "NeighborSearcher.h"
#include <map>

namespace MeshFunc{
void Mesh_To_Level_Set(const SegmentMesh<2>& mesh, const Grid<2>& grid, LevelSet<2>& level_set) {
	Array<Vector2> boundary_ps;
	std::map<int, int> boundary_ps_to_seg_index;	//the mapping between the index of boundary pos and segment element index

	const Array<Vector2i>& es= mesh.Elements();
	const Array<Vector2>& vs=mesh.Vertices();

	//sample points on segments, include p0 but not p1
	//make sure that every grid cell on the boundary has one point sampled
	for (int i=0;i<es.size();i++) {
		int p0_index = es[i][0];
		int p1_index = es[i][1];
		Vector2 p0_pos = vs[p0_index];
		Vector2 p1_pos = vs[p1_index];
		Vector2 v01_hat=(p1_pos-p0_pos).normalized();
		real step_size=std::min(abs(v01_hat.dot(Vector2::Unit(0)*grid.dx)),abs(v01_hat.dot(Vector2::Unit(1)*grid.dx)));
		Vector2 step_unit =step_size*v01_hat;
		int num_step = (int)((p1_pos-p0_pos).norm()/step_size);

		for (int j = 0; j < num_step; j++) {
			boundary_ps_to_seg_index[(int)boundary_ps.size()]=i;
			boundary_ps.push_back(p0_pos+j*step_unit);}}

	NeighborKDTree<2> nbs_searcher = NeighborKDTree<2>();
	nbs_searcher.Build_Data(boundary_ps);
	// set phi for level set
	int cell_num=grid.cell_counts.prod();
	for(int i=0;i<cell_num;i++){
		const Vector2i cell=grid.Cell_Coord(i);
		Vector2 pos=grid.Center(cell);
		int nb=nbs_searcher.Find_Nearest_Nb(grid.Center(cell));
		int seg_index = boundary_ps_to_seg_index[nb];
		Segment<2> seg(vs[es[seg_index][0]],vs[es[seg_index][1]]);
		Vector2 vec = seg.p1-seg.p0;
		Vector2 normal = Vector2(-vec[1],vec[0]); //rotate the segment vec 90 degrees counter clockwise
		bool is_negative=signbit(normal.dot(Segment<2>::Closest_Point_Vector(seg.p0,seg.p1,pos)));
		level_set.phi(cell)=is_negative?seg.Phi(pos):-seg.Phi(pos);}
}

void Mesh_To_Level_Set(TriangleMesh<3>& mesh, const Grid<3>& grid, LevelSet<3>& level_set) {
	real epsilon = grid.dx;
	int cell_num=grid.cell_counts.prod();
	Update_Normals(mesh);
	Array<Array<int>> neighbor_hashing(cell_num);
	const Array<Vector3i>& es= mesh.Elements();
	const Array<Vector3>& vs=mesh.Vertices();
	const Array<Vector3>& vn=*mesh.Normals();
	Field<int,3> sign_field(grid.cell_counts,-10);

	//Store triangle index to the neighboring cell for neighbor_hashing
	for (int i=0;i<es.size();i++) {
		int p0_index = es[i][0];
		int p1_index = es[i][1];
		int p2_index = es[i][2];
		Vector3 p0_pos = vs[p0_index];
		Vector3 p1_pos = vs[p1_index];
		Vector3 p2_pos = vs[p2_index];
		Array<Vector3> ps={p0_pos,p1_pos,p2_pos};
		Box<3> bounding_box=Bounding_Box<3>(ps);
		Box<3> enlarged_box=bounding_box.Enlarged(Vector3::Ones()*epsilon);
		Vector3i min_coord = grid.Clamped_Cell_Coord(grid.Cell_Coord(enlarged_box.min_corner));
		Vector3i max_coord = grid.Clamped_Cell_Coord(grid.Cell_Coord(enlarged_box.max_corner));
		for (int x = min_coord[0]; x <= max_coord[0]; x++) {
			for (int y = min_coord[1]; y <= max_coord[1]; y++) {
				for (int z = min_coord[2]; z <= max_coord[2]; z++) {
					int cell_index=grid.Cell_Index(Vector3i(x,y,z));
					neighbor_hashing[cell_index].push_back(i);}}}
	}
	
	// find the cloest triangle from neighbor hashing
	for (int i = 0; i < cell_num; i++) {
		Array<int>& neighbors=neighbor_hashing[i];
		Vector3i cell=grid.Cell_Coord(i);
		Vector3 pos=grid.Center(cell);
		if(neighbors.size()==0){continue;}
		real min_dist=FLT_MAX; int min_triangle_index; bool is_negative;
		for (int j = 0; j < neighbors.size(); j++) {
			int triangle_index=neighbors[j];
			Triangle<3> triangle(vs[es[triangle_index][0]],vs[es[triangle_index][1]],vs[es[triangle_index][2]]);
			Triangle<3>::ClosestPointFeature feature;
			Vector3 closest_point_vector=Triangle<3>::Closest_Point_Vector(triangle.p0,triangle.p1,triangle.p2,pos,feature);
			real dist=closest_point_vector.norm();
			if(dist<min_dist){
				Vector3 normal;
				switch (feature.m_Type) {
				case(Triangle<3>::ClosestPointFeature::Type::Face): { //use face normal
					normal=Normal(triangle.p0,triangle.p1,triangle.p2);is_negative=signbit(normal.dot(closest_point_vector));
				}break;
				case(Triangle<3>::ClosestPointFeature::Type::Edge): { //use averaged vertex normal
					int edge=feature.m_Index;
					int vertex_start=Triangle<3>::ClosestPointFeature::c_EdgeStartPoint[edge];
					int vertex_end=Triangle<3>::ClosestPointFeature::c_EdgeEndPoint[edge];
					normal=(real)0.5*(vn[es[triangle_index][vertex_start]]+vn[es[triangle_index][vertex_end]]);
				}break;
				case(Triangle<3>::ClosestPointFeature::Type::Vertex): { //use vertex normal
					normal=vn[es[triangle_index][feature.m_Index]];is_negative=signbit(normal.dot(closest_point_vector));
				}break;}
				is_negative=signbit(normal.dot(closest_point_vector));	//calculate the sign
				min_dist=dist;min_triangle_index=triangle_index;}}
		
		level_set.phi(cell)=is_negative?-min_dist:min_dist;
		sign_field(cell)=is_negative?-1:1;}

	//Fill the sign of phi for fast marching
	Flood_Fill(sign_field,grid,0,-1,1);
	for (int i = 0; i < cell_num; i++) {Vector3i cell=grid.Cell_Coord(i);level_set.phi(cell)=sign_field(cell)*abs(level_set.phi(cell));}
	
	//Use fast marching to calculate the area without a real phi value
	level_set.Fast_Marching();
}
};
#endif
