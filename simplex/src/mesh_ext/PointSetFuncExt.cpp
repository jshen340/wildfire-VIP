#include "PointSetFuncExt.h"

namespace PointSetFunc{
        real Initialize_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles, const real randomness)
    {
        TriangleMesh<3> mesh;
        MeshFunc::Initialize_Sphere_Mesh(R, &mesh, sub);
        MeshFunc::Update_Normals(mesh);
        particles.Resize((int)mesh.Vertices().size());
        for (int i = 0; i < particles.Size(); i++) {
            particles.X(i) = mesh.Vertices()[i] + c;
            particles.X(i)[0] += randomness * 2 * (((real)rand() / (RAND_MAX + 1.)) - 0.5);
            particles.X(i)[1] += randomness * 2 * (((real)rand() / (RAND_MAX + 1.)) -0.5);
            particles.X(i)[2] += randomness * 2 * (((real)rand() / (RAND_MAX + 1.)) -0.5);
            particles.X(i) = particles.X(i).normalized();
            particles.X(i) *= R;
            particles.V(i) = Vector3::Zero();
            particles.F(i) = Vector3::Zero();
            particles.M(i) = (real)1;
            Vector3 normal = (*mesh.Normals())[i];
            Vector3 t1 = -Orthogonal_Vector(normal);
            Vector3 t2 = t1.cross(normal);
            particles.E(i).col(0) = t1;
            particles.E(i).col(1) = t2;
            particles.E(i).col(2) = normal;
        }
        real dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
        return dx;
    }

    real Initialize_Half_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles, std::vector<int>* is_boundary)
{
	if (is_boundary!=nullptr) is_boundary->clear();
	TriangleMesh<3> mesh;
	MeshFunc::Initialize_Sphere_Mesh(R, &mesh, sub);
	MeshFunc::Update_Normals(mesh);
	int size = 0;
	real cutoff = -0.3 * R;
	real scale = R / (sqrt(R * R - cutoff * cutoff));
	for (int i = 0; i < mesh.Vertices().size(); i++) {
		if (mesh.Vertices()[i][1] < cutoff) continue;
		else size++;
	}
	particles.Resize(size);
	int j = 0;
	for (int i = 0; i < mesh.Vertices().size(); i++) {
		particles.X(j) = mesh.Vertices()[i] + c;
		particles.X(j)[0] += R / 20. * (((real)rand() / (RAND_MAX + 1.)) - 0.5);
		particles.X(j)[1] += R / 20. * (((real)rand() / (RAND_MAX + 1.)) - 0.5);
		particles.X(j)[2] += R / 20. * (((real)rand() / (RAND_MAX + 1.)) - 0.5);
		particles.X(j) = particles.X(j).normalized();
		particles.X(j) *= R;
		if (particles.X(j)[1] < cutoff) continue;
		particles.X(j)[1] -= cutoff;
		particles.X(j)[0] *= scale;
		particles.X(j)[2] *= scale;
		if (particles.X(j)[1] < 0.) continue;
		particles.V(j) = Vector3::Zero();
		particles.F(j) = Vector3::Zero();
		particles.M(j) = (real)1;
		if (is_boundary != nullptr) {
			int b = 0;
			b = particles.X(j)[1] < 1e-8 ? 1 : 0;
			is_boundary->push_back(b);
		}
		Vector3 normal = (*mesh.Normals())[i];
		Vector3 t1 = -Orthogonal_Vector(normal);
		Vector3 t2 = t1.cross(normal);
		particles.E(j).col(0) = t1;
		particles.E(j).col(1) = t2;
		particles.E(j).col(2) = normal;
		j++;
	}
	particles.Resize(j);
	real dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());

	std::vector<Vector3, Eigen::aligned_allocator<Vector3>> pts;

	int boundary_multiplier = 1;
	int num_boundary = int(boundary_multiplier * 2 * pi * R / dx);
	for (int i = 0; i < num_boundary; i++) {
		real x_coord = R * cos(((real)i / num_boundary) * (2 * pi));
		real z_coord = R * sin(((real)i / num_boundary) * (2 * pi));
		Vector3 curr_x;
		curr_x << x_coord, 0, z_coord;
		pts.push_back(curr_x);
	}

	int prev_size = particles.Size();
	particles.Resize(particles.Size() + (int)pts.size());
	for (int i = 0; i < pts.size(); i++) {
		particles.X(prev_size + i) = pts[i] + c;
		particles.V(prev_size + i) = Vector3::Zero();
		particles.F(prev_size + i) = Vector3::Zero();
		particles.M(prev_size + i) = (real)1;
		Vector3 normal = (particles.X(prev_size + i) - c).normalized();
		Vector3 t1 = -Orthogonal_Vector(normal);
		Vector3 t2 = t1.cross(normal);
		particles.E(prev_size + i).col(0) = t1;
		particles.E(prev_size + i).col(1) = t2;
		particles.E(prev_size + i).col(2) = normal;
		is_boundary->push_back(1);
	}
	return dx;
}

real Add_Sphere_Points(const Vector3& c, const real R, const int sub, GeometryParticles<3>& particles)
{
	TriangleMesh<3> mesh;
	MeshFunc::Initialize_Sphere_Mesh(R, &mesh, sub);
	MeshFunc::Update_Normals(mesh);
	int osize=particles.Size();
	particles.Resize(osize+(int)mesh.Vertices().size());
	for (int i = osize; i < particles.Size(); i++) {
		particles.X(i) = mesh.Vertices()[i-osize] + c;
		particles.V(i) = Vector3::Zero();
		particles.F(i) = Vector3::Zero();
		particles.M(i) = (real)1;
		Vector3 normal = (*mesh.Normals())[i-osize];
		Vector3 t1 = -Orthogonal_Vector(normal);
		Vector3 t2 = t1.cross(normal);
		particles.E(i).col(0) = t1;
		particles.E(i).col(1) = t2;
		particles.E(i).col(2) = normal;
	}
	real dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
	return dx;
}


real Initialize_Circle_Points(const Vector3& c, const real R, const Vector3& normal, const real dx, GeometryParticles<3>& particles)
{
	int p_num = (int)ceil(two_pi * R / dx);
	TriangleMesh<3> mesh;
	MeshFunc::Initialize_Circle_Mesh(c, R, normal, &mesh, p_num);
	MeshFunc::Update_Normals(mesh);
	particles.Resize((int)mesh.Vertices().size());
	for (int i = 0; i < particles.Size(); i++) {
		particles.X(i) = mesh.Vertices()[i];
		particles.V(i) = Vector3::Zero();
		particles.F(i) = Vector3::Zero();
		particles.M(i) = (real)1;
		Vector3 normal = (*mesh.Normals())[i];
		Vector3 t1 = -Orthogonal_Vector(normal);
		Vector3 t2 = t1.cross(normal);
		particles.E(i).col(0) = t1;
		particles.E(i).col(1) = t2;
		particles.E(i).col(2) = normal;
	}
	real avg_dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
	return avg_dx;
}

real Initialize_Ellipsoid_Points(const Vector3& C, const real R, const int sub, const real a, const real b, const real c,  GeometryParticles<3>& particles)
{
	TriangleMesh<3> mesh;
	MeshFunc::Initialize_Ellipsoid_Mesh(R, &mesh, a,b,c, sub);
	MeshFunc::Update_Normals(mesh);
	particles.Resize((int)mesh.Vertices().size());
	for (int i = 0; i < particles.Size(); i++) {
		particles.X(i) = mesh.Vertices()[i] + C;
		particles.V(i) = Vector3::Zero();
		particles.F(i) = Vector3::Zero();
		particles.M(i) = (real)1;
		Vector3 normal = (*mesh.Normals())[i];
		Vector3 t1 = -Orthogonal_Vector(normal);
		Vector3 t2 = t1.cross(normal);
		particles.E(i).col(0) = t1;
		particles.E(i).col(1) = t2;
		particles.E(i).col(2) = normal;
	}
	real dx = MeshFunc::Average_Edge_Length<3>(mesh.Vertices(), mesh.Elements());
	return dx;
}


}