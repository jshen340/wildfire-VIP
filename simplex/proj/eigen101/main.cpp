//#####################################################################
// Eigen 101
// Dartmouth CS89.18/189
//#####################################################################
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

////Aliasing class names 
using T=double;	////Attn: using (c++11 feature, equivalent to typedef)
template<class TYPE,int d> using Vector=Eigen::Matrix<TYPE,d,1>;
template<class TYPE> using VectorN=Eigen::Matrix<TYPE,-1,1>;	////TODO: why -1?
using Vector2=Vector<T,2>;
using Vector3=Vector<T,3>;
using Vector4=Vector<T,4>;
using Vector3i=Vector<int,3>;
using Matrix2=Eigen::Matrix<T,2,2>;
using Matrix3=Eigen::Matrix<T,3,3>;
using VectorX=Eigen::Matrix<T,Eigen::Dynamic,1>;
using MatrixX=Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template<class TYPE> using SparseMatrix=Eigen::SparseMatrix<TYPE,Eigen::RowMajor,int>;
template<class TYPE> using Triplet=Eigen::Triplet<TYPE,int>;
using CG=Eigen::ConjugateGradient<SparseMatrix<T>,Eigen::Upper,Eigen::IdentityPreconditioner>;
using InnerIterator=SparseMatrix<T>::InnerIterator;

void Test_Name(std::string name)
{
	std::cout<<"---------- "<<name<<" ----------"<<std::endl;
}

void Dense_Tests()
{
	Test_Name("Initialization");
	Vector2 v2(1,2);
    Vector3 v3(1,2,3);
	Vector4 v4;v4<<1,2,3,4;
	using Vector8=Vector<T,8>;	////define your own type
	Vector8 v8=Vector8::Ones();
	Matrix3 m3=Matrix3::Ones();
	m3<<1,2,3,4,5,6,7,8,9;	////by default a matrix storage is column-major
	
	std::cout<<"v2= "<<v2.transpose()
		<<"\nv3= "<<v3.transpose()
		<<"\nv4= "<<v4.transpose()
		<<"\nv8= "<<v8.transpose()
		<<"\nm3:\n"<<m3<<std::endl;

	Test_Name("Type conversion");
    Vector3 vd(1.,2.,3.);
    Vector3i vi=vd.cast<int>();
    Vector3 vt=vi.cast<T>();
    
	std::cout<<"vi= "<<vi.transpose()
		<<"\nvt= "<<vt.transpose()<<std::endl;

	Test_Name("Dynamic size");
	int r=3;int c=4;
	MatrixX mx(r,c);
	int ct=0;
	for(auto i=0;i<mx.rows();i++){
		for(auto j=0;j<mx.cols();j++){
			mx(i,j)=(T)(++ct);}}
	int n=4;
	VectorX vx(n);
	for(auto i=0;i<vx.rows();i++)vx[i]=(T)(i+1);
	
	std::cout<<"mx:\n"<<"rows= "<<mx.rows()<<", cols= "<<mx.cols()<<"\n"<<mx<<std::endl;
	std::cout<<"\nvx:\nsize= "<<vx.rows()<<"\n"<<vx.transpose()<<std::endl;

	Test_Name("Arithmetics");
	VectorX vx2=mx*vx;	////matrix-vector multiplication
	Matrix3 i3=Matrix3::Identity();
	T q=v3.transpose()*i3*v3;	////vector-matrix-vector multiplication
	MatrixX mx2=mx.transpose()*i3;	////matrix-matrix multiplication
	T dot=vx.dot(vx2);	////dot product
	VectorX cp=vx.cwiseProduct(vx2);
	T sum=vx.sum();		////reduction operations: sum, prod, min, max, norms
	T prod=vx.prod();
	T min_coef=vx.minCoeff();
	T max_coef=vx.maxCoeff();
	T norm=vx.norm();
	T sq_norm=vx.squaredNorm();
	T L1_norm=vx.lpNorm<1>();
	T L2_norm=vx.lpNorm<2>();
	T Linf_norm=vx.lpNorm<Eigen::Infinity>();

	std::cout<<"mx*vx= "<<vx2.transpose()
		<<"\nvTMv= "<<q
		<<"\nmx*i3=\n"<<mx2
		<<"\nvx dot vx2= "<<dot
		<<"\ncomponent wise product= "<<cp.transpose()
		<<"\nsum= "<<sum
		<<"\nprod= "<<prod
		<<"\nmin_coef= "<<min_coef
		<<"\nmax_coef= "<<max_coef
		<<"\nnorm= "<<norm
		<<"\nsquared norm= "<<sq_norm
		<<"\nL1 norm= "<<L1_norm
		<<"\nL2 norm= "<<L2_norm
		<<"\nLinf norm= "<<Linf_norm<<std::endl;

	Test_Name("Block");
	Matrix2 b=mx.block<2,2>(0,1);	////a block with starting index (0,1) and size (2,2)
	Matrix3 b2=mx.block(0,1,3,2);	////a block with starting index (0,1) and size (3,2)
	VectorX br=mx.row(0);
	VectorX bc=mx.col(0);
	std::cout<<"block2x2=\n"<<b<<std::endl;
	std::cout<<"block3x2=\n"<<b2<<std::endl;
	std::cout<<"row= "<<br.transpose()<<std::endl;
	std::cout<<"col= "<<bc.transpose()<<std::endl;

	Test_Name("Map");
	std::vector<T> array={1,2,3,4,5,6,7,8};
	Eigen::Map<Vector2> mp(&array[1]);
	Eigen::Map<Eigen::Matrix<T,2,4,Eigen::RowMajor> > mp2(&array[0]);
	std::cout<<"mp= "<<mp.transpose()<<std::endl;
	std::cout<<"mp2=\n"<<mp2<<std::endl;
}

void Sparse_Tests()
{
	Test_Name("Construct sparse matrix");
	int nc=3;
	SparseMatrix<T> A(nc,nc);
	std::vector<Triplet<T> > elements={{0,0,2.},{0,1,-1.},{1,0,-1.},{1,1,2.},{1,2,-1.},{2,2,2.},{2,1,-1.}};
	A.setFromTriplets(elements.begin(),elements.end());
	A.makeCompressed();
	std::cout<<"A=\n"<<A;

	Test_Name("Element Access");
	T a11=A.coeff(1,1);	////read only
	std::cout<<"A[1,1]= "<<a11<<std::endl;
	A.coeffRef(1,1)=3;	////modify value
	std::cout<<"A[1,1]= "<<A.coeff(1,1)<<std::endl;

	T a02=A.coeff(0,2);
	std::cout<<"A[0,2]= "<<a02<<std::endl;
	A.coeffRef(0,2)=4;	////modify zero value to nonzero, expensive
	std::cout<<"A[0,2]= "<<A.coeff(0,2)<<std::endl;

	Test_Name("Solve sparse system");
	VectorN<T> x(nc);x.fill((T)0);
	VectorN<T> b(nc);b<<1.,0.,1.;
	CG cg;
	cg.compute(A);
	x=cg.solve(b);
	std::cout<<"x= "<<x.transpose()<<std::endl;
	std::cout<<"b= "<<b.transpose()<<std::endl;

	Test_Name("Interate nonzero elements");
	std::cout<<"nnz= "<<A.nonZeros()<<std::endl;
	for(int r=0;r<A.rows();r++){
		for(InnerIterator iter(A,r);iter;++iter){
			std::cout<<"A["<<iter.row()<<","<<iter.col()<<"]= "<<iter.value()<<std::endl;}}
	
	Test_Name("Raw pointers");
	T* val=A.valuePtr();
	int* col=A.innerIndexPtr();
	int* ptr=A.outerIndexPtr();
	std::cout<<"val= ";
	for(int i=0;i<A.nonZeros();i++){
		std::cout<<val[i]<<" ";}
	std::cout<<"\ncol= ";
	for(int i=0;i<A.nonZeros();i++){
		std::cout<<col[i]<<" ";}
	std::cout<<"\nptr= ";
	for(int i=0;i<A.rows()+1;i++){
		std::cout<<ptr[i]<<" ";}
}

void Sparse_Test2()
{
	SparseMatrix<T> A;
	int r=3;int c=3;
	A.resize(r,c);
	A.resizeNonZeros(7);
	std::cout<<"A: "<<A.rows()<<", "<<A.cols()<<": "
		<<A.data().size()<<", "<<A.nonZeros()<<", "
		<<A.innerSize()<<", "<<A.outerSize()<<std::endl;
	T* val=A.valuePtr();
	int* col=A.innerIndexPtr();
	int* ptr=A.outerIndexPtr();
	
	std::vector<int> ptr_array={0,2,5,7};
	std::vector<int> col_array={0,1,0,1,2,1,2};
	std::vector<T> val_array={2.,-1.,-1.,2.,-1.,-1.,2.};
	for(auto i=0;i<ptr_array.size();i++)ptr[i]=ptr_array[i];
	for(auto i=0;i<col_array.size();i++){col[i]=col_array[i];val[i]=val_array[i];}
	std::cout<<"A:\n"<<A<<std::endl;
}

void Transformation_Tests()
{
	////TODO
}

int main()
{
	//Dense_Tests();
	//Sparse_Tests();
	//Transformation_Tests();
	Sparse_Test2();
    return 0;
}