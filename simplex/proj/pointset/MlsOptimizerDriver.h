#ifndef __MlsOptimizerDriver_h__
#define __MlsOptimizerDriver_h__
#include "Driver.h"
#include "LeastSquares.h"
#include "Mesh.h"

#include "PointSet.h"
#include "PointSetFunc.h"

// #include <autodiff/forward.hpp>
// #include <autodiff/forward/eigen.hpp>
// #include <autodiff/reverse.hpp>
// #include <autodiff/reverse/eigen.hpp>
// using namespace autodiff;

template<int d> class MlsOptimizerDriver : public Driver
{public:
	using Base=Driver;
	using VectorDi=Vector<int,d>;using VectorD=Vector<real,d>;
	using VectorTi=Vector<int,d-1>;using VectorT=Vector<real,d-1>;
    using MatrixD=Matrix<real,d>;

	PointSet<d> ps;
	Array<real> data,approx;
	Vector2 domain=Vector2(-1.,1.);
	int count=16;

	virtual void Initialize()
	{
		using namespace LeastSquares;
		switch(test){
		case 1: {////MLS<2,2>
			if constexpr (d==3) {
                real dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
                ps.Initialize(dx,3);
                ps.Update();
                //// move to zero
                VectorD pv=ps.points->X(0);
                for(int i=0;i<ps.points->Size();i++)
                    ps.points->X(i)-=pv;
                ps.Update();
			}
		}break;
        case 2: {////MLS<2,2>
			if constexpr (d==3) {
                real dx=PointSetFunc::Initialize_Sphere_Points(VectorD::Zero(),(real)1,3,*ps.points);
                PointSetFunc::Add_Sphere_Points(VectorD(1,1,2),(real)1,3,*ps.points);
                std::cout<<"PSize"<<ps.points->Size()<<std::endl;
                ps.Initialize(dx,3);
                ps.Update();
                ps.Reinitialize_Local_Frames();
                // for(int i=0;i<ps.points->Size();i++){
                //     std::cout<<ps.points->X(i).transpose()<<std::endl;
                // }
                Ray_Cast_Pointset();

			}
		}break;
		default:{////WLS<2,1>
		}break;
		}
	}

    virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		switch(test){
		case 1:{
		    Optimize_Pointset();
		}break;
        case 2:{
		}break;
        default:{////WLS<2,1>
		}break;
		}
        
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Base::Write_Output_Files(frame);

		{std::string file_name=frame_dir+"/tracker_circles";
		PointSetFunc::Write_Tracker_Circles_To_File<d>(file_name,*ps.points);}

		{std::string file_name=frame_dir+"/segment_mesh";
		PointSetFunc::Write_Local_Frames_To_File(file_name,*ps.points);}
		
		{std::string file_name = frame_dir + "/point_force";
		Write_Segments_To_File_3d_Fast<d, real>(ps.points->XRef(), ps.points->FRef(), file_name); }
		{std::string file_name = frame_dir + "/point_velocity";
		Write_Segments_To_File_3d_Fast<d, real>(ps.points->XRef(), ps.points->VRef(), file_name); }

		std::cout<<"Write to frame "<<frame<<std::endl;
	}

	void Optimize_Pointset()
    {
        ps.Update();
        ps.Reinitialize_Local_Frames();
        const int d=3;
        int p=0; //// optimize at the implicit surface on p
        
        int n=ps.tang_nbs[p].size();
        VectorX data(n*d);
        VectorX weights(n);
		VectorD x_p=ps.points->X(p);
        MatrixD local_frame=ps.points->E(p);
        std::function<real(const int)> Point_Position_Z=[&](const int idx)->real
				{return (ps.points->X(idx).transpose()*local_frame)[2];};

        //// print surface coefficient
		VectorX coeff=ps.Fit_Local_MLS(0,Point_Position_Z);
        std::cout<<"Surface Coefficient: "<<coeff.transpose()<<std::endl;

        //// neighboring particle data in local coordinate
		int pi=0;
        for(int i=0;i<n;i++){
			int nb=ps.tang_nbs[p][i];
			if(nb==p) {pi=i; weights(i)=(real)1;}
			else weights(i)=(real)1/(real)n;
            VectorD th;ps.Project_To_TPlane(ps.points->X(nb)-x_p,local_frame,th);
			for(int j=0;j<d;j++)data(i*d+j)=th[j];}
        
        //// compute gradient of the curvature
        for(int i=0;i<n;i++){
            int nb=ps.tang_nbs[p][i];
            ps.points->V(nb)=VectorD::Zero();
            std::pair<MatrixX,VectorD> pair=Compute_Gradient(data, weights, pi, i);
            
            //// grad on coefficients
            MatrixX grad_coef=pair.first;
            //// par_(curv)=c3c3)/par_x
            VectorD local_grad=coeff(3)*grad_coef.col(3)+coeff(4)*grad_coef.col(4);
            ps.points->V(nb)-=local_frame*local_grad;

            //// grad on MLS objective
            VectorD grad_mlsobj=pair.second;
            ps.points->V(nb)-=local_frame*grad_mlsobj*100.0f;
        }

        //// update velocity
        for(int i=0;i<n;i++){
            int nb=ps.tang_nbs[p][i];
            ps.points->X(nb)+=ps.points->V(nb)*0.02;}

        ps.Update();
    }

    //// compute partial beta/partial x manualy
    std::pair<MatrixX,VectorD> Compute_Gradient(const VectorX& data, const VectorX& weights, const int p, const int j)
	{
        ///// assumption: par w/par x=0;
        const int d=3;
        MatrixX B(weights.size(),6);
        MatrixX W(weights.size(),weights.size());
		VectorX f(weights.size());

        //// set matrix B and W
		for(int i=0;i<weights.size();i++){
            real w=sqrt(weights(i));
            W.coeffRef(i,i)=weights(i);
            B.coeffRef(i,0)=(real)1;
            B.coeffRef(i,1)=data[i*d];
            B.coeffRef(i,2)=data[i*d+1];
            B.coeffRef(i,3)=pow(data[i*d],2);
            B.coeffRef(i,4)=pow(data[i*d+1],2);
            B.coeffRef(i,5)=data[i*d]*data[i*d+1];
			f[i]=data[i*d+2];
        }

        //// set partial derivative
        MatrixX parB_parx[d-1];
        for(int i=0;i<d-1;i++){parB_parx[i].resize(weights.size(),6);parB_parx[i].fill(0);}
        
        parB_parx[0].coeffRef(j,1)=1;
        parB_parx[0].coeffRef(j,3)=2.0*data[j*3];
        parB_parx[0].coeffRef(j,5)=data[j*3+1];

        parB_parx[1].coeffRef(j,2)=1;
        parB_parx[1].coeffRef(j,4)=2.0*data[j*3+1];
        parB_parx[1].coeffRef(j,5)=data[j*3];


        MatrixX BtWB=B.transpose()*W*B;
        MatrixX BtWB_inv=BtWB.inverse();

        //// par(BtWB)/par(x)=par(Bt)/par(x)*W*B + Bt*W*par(B)/par(x)
        MatrixX parBtWB_parx0=parB_parx[0].transpose()*W*B+B.transpose()*W*parB_parx[0];
        MatrixX parBtWB_parx1=parB_parx[1].transpose()*W*B+B.transpose()*W*parB_parx[1];

        //// gradient on coefficient
        MatrixX grad_coef(d,6);
        //// par((BtWB)^(-1)*BtWz)/par(x)=-(BtWB)^(-1)*par(BtWB)/par(x)*(BtWB)^(-1)*BtWz+(BtWB)^(-1)*par(BtWz)/par(x)
        grad_coef.row(0)=(-BtWB_inv*parBtWB_parx0*BtWB_inv*B.transpose()*W*f+BtWB_inv*parB_parx[0].transpose()*W*f).transpose();
        //// par((BtWB)^(-1)*BtWz)/par(y)=-(BtWB)^(-1)*par(BtWB)/par(y)*(BtWB)^(-1)*BtWz+(BtWB)^(-1)*par(BtWz)/par(y)
        grad_coef.row(1)=(-BtWB_inv*parBtWB_parx1*BtWB_inv*B.transpose()*W*f+BtWB_inv*parB_parx[1].transpose()*W*f).transpose();
        //// par((BtWB)^(-1)*BtWz)/par(z)=(BtWB)^(-1)*par(BtWz)/par(z)
        grad_coef.row(2)=(BtWB_inv*B.transpose()*W).col(j).transpose();
        // return grad;

        //// gradient on MLS objective
        VectorX coeff=(W*B).colPivHouseholderQr().solve(W*f);
        VectorD grad_mlsobj;
        grad_mlsobj[0]=(2*W*(B*coeff-f)).dot(parB_parx[0]*coeff);
        grad_mlsobj[1]=(2*W*(B*coeff-f)).dot(parB_parx[1]*coeff);
        grad_mlsobj[2]=-(2*W*(B*coeff-f))(j);

        return std::pair<MatrixX,VectorD>(grad_coef,grad_mlsobj);
	}

	void Update_Data(int n, int deg)
	{
		data.clear();
		for(int i=0;i<n;i++){
			for(int j=0;j<=deg;j++){
				data.push_back(particles.X(i)[j]);}}
	}



	void Ray_Cast_Pointset()
    {
        Array<VectorD> rays;
        Array<real> depths;
        Array<MatrixD> lfs;
        Array<VectorX> cs;
        VectorD eye(0,0,-3.5);
        Get_Ray(64, eye, VectorD(0,0,0), VectorD(0,1,0), 2.2, rays);
        depths.resize(rays.size());
        lfs.resize(rays.size());
        cs.resize(rays.size());
        
        //// Compute ray cast
        for(int i=0;i<rays.size();i++){
            depths[i]=ps.Ray_Casting_MLS(eye, rays[i], 10, lfs[i], cs[i]);
            VectorD p=eye+depths[i]*rays[i];
            std::cout<<"["<<i/64<<","<<i%64<<"]"<<"p=("<<p.transpose()<<")"<<p.norm()<<" ray=("<<rays[i].transpose()<<"), depth="<<depths[i]<<std::endl;
            std::cout<<"["<<i/64<<","<<i%64<<"]"<<"("<<lfs[i]<<")"<<std::endl;
        }

        //// Save ray-casting depth in binary format
        {
            std::ofstream ofs("./depth_array.dat",std::ios::out|std::ios::binary);
            for(int i=0;i<rays.size();i++)
                ofs.write((char*)(&depths[i]), sizeof(real));
            ofs.close();
        }
        {
            std::ofstream ofs("./coef_array.dat",std::ios::out|std::ios::binary);
            
            for(int i=0;i<cs.size();i++){
                VectorX coef=cs[i];
                ofs.write((char*)(&coef(0)), sizeof(real));
                ofs.write((char*)(&coef(1)), sizeof(real));
                ofs.write((char*)(&coef(2)), sizeof(real));
                ofs.write((char*)(&coef(3)), sizeof(real));
                ofs.write((char*)(&coef(4)), sizeof(real));
                ofs.write((char*)(&coef(5)), sizeof(real));}
            ofs.close();
        }
        {
            std::ofstream ofs("./ray_array.dat",std::ios::out|std::ios::binary);
            if constexpr(d==2){
                for(int i=0;i<rays.size();i++){
                    VectorD ray=rays[i];
                    ofs.write((char*)(&ray(0)), sizeof(real));
                    ofs.write((char*)(&ray(1)), sizeof(real));}}
            else if constexpr(d==3){
                for(int i=0;i<rays.size();i++){
                    VectorD ray=rays[i];
                    ofs.write((char*)(&ray(0)), sizeof(real));
                    ofs.write((char*)(&ray(1)), sizeof(real));
                    ofs.write((char*)(&ray(2)), sizeof(real));}}
            ofs.close();
        }
        {
            std::ofstream ofs("./hit_array.dat",std::ios::out|std::ios::binary);
            if constexpr(d==2){
                for(int i=0;i<rays.size();i++){
                    VectorD p=eye+depths[i]*rays[i];
                    ofs.write((char*)(&p(0)), sizeof(real));
                    ofs.write((char*)(&p(1)), sizeof(real));}}
            else if constexpr(d==3){
                for(int i=0;i<rays.size();i++){
                    VectorD p=eye+depths[i]*rays[i];
                    ofs.write((char*)(&p(0)), sizeof(real));
                    ofs.write((char*)(&p(1)), sizeof(real));
                    ofs.write((char*)(&p(2)), sizeof(real));}}
            ofs.close();
        }
        {
            std::ofstream ofs("./lf_array.dat",std::ios::out|std::ios::binary);
            if constexpr(d==2){
                for(int i=0;i<rays.size();i++){
                    MatrixD lf=lfs[i];
                    ofs.write((char*)(&lf(0,0)), sizeof(real));
                    ofs.write((char*)(&lf(0,1)), sizeof(real));
                    ofs.write((char*)(&lf(1,0)), sizeof(real));
                    ofs.write((char*)(&lf(1,1)), sizeof(real));}}
            else if constexpr(d==3){
                for(int i=0;i<rays.size();i++){
                    MatrixD lf=lfs[i];
                    ofs.write((char*)(&lf(0,0)), sizeof(real));
                    ofs.write((char*)(&lf(0,1)), sizeof(real));
                    ofs.write((char*)(&lf(0,2)), sizeof(real));
                    ofs.write((char*)(&lf(1,0)), sizeof(real));
                    ofs.write((char*)(&lf(1,1)), sizeof(real));
                    ofs.write((char*)(&lf(1,2)), sizeof(real));
                    ofs.write((char*)(&lf(2,0)), sizeof(real));
                    ofs.write((char*)(&lf(2,1)), sizeof(real));
                    ofs.write((char*)(&lf(2,2)), sizeof(real));}}
            ofs.close();
        }
        {
            std::ofstream ofs("./point_array.dat",std::ios::out|std::ios::binary);
            if constexpr(d==2){
                for(int i=0;i<ps.points->Size();i++){
                    const VectorD p=ps.points->X(i);
                    ofs.write((char*)(&p(0)), sizeof(real));
                    ofs.write((char*)(&p(1)), sizeof(real));}}
            else if constexpr(d==3){
                for(int i=0;i<ps.points->Size();i++){
                    const VectorD p=ps.points->X(i);
                    ofs.write((char*)(&p(0)), sizeof(real));
                    ofs.write((char*)(&p(1)), sizeof(real));
                    ofs.write((char*)(&p(2)), sizeof(real));}}
            ofs.close();
        }
        {
            std::ofstream ofs("./norm_array.dat",std::ios::out|std::ios::binary);
            if constexpr(d==2){
                for(int i=0;i<rays.size();i++){
                    VectorD p=eye+depths[i]*rays[i]; VectorD n=ps.Normal(p);
                    ofs.write((char*)(&n(0)), sizeof(real));
                    ofs.write((char*)(&n(1)), sizeof(real));}}
            else if constexpr(d==3){
                for(int i=0;i<rays.size();i++){
                    VectorD p=eye+depths[i]*rays[i]; VectorD n=ps.Normal(p);
                    ofs.write((char*)(&n(0)), sizeof(real));
                    ofs.write((char*)(&n(1)), sizeof(real));
                    ofs.write((char*)(&n(2)), sizeof(real));}}
            ofs.close();
        }
        {
            std::ofstream ofs("./r.dat",std::ios::out|std::ios::binary);
            ofs.write((char*)(&ps.v_r), sizeof(real));
            ofs.close();
        }
    }

    void Get_Ray(const int n, const VectorD eye, const VectorD look, const VectorD up, const real d, Array<VectorD>& rays){
        real step = (real)1./(n-1);

        VectorD w=(look-eye).normalized();
        VectorD u=(w.cross(up)).normalized();
        VectorD v=(u.cross(w)).normalized();

        rays.resize(n*n);
        for(int yi=0;yi<n;yi++)
        for(int xi=0;xi<n;xi++){
            real x=xi*step*2-1;
            real y=yi*step*2-1;
            rays[yi*n+xi]=(x*u+y*v+d*w).normalized();
        }
        return;
    }
};

#endif