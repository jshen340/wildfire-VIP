#include "io_adapter.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{
	using namespace IOAdapter;
	std::string output_dir="output";

	////grid properties
	int cell_counts[3]={16,16,16};
	double dx=1./(double)cell_counts[0];
	double domain_min[3]={0.,0.,0.};

	////particle property
	double px[6]={0.,0.,0., .5,.5,.5};
	int pn=2;

	int vn=cell_counts[0]*cell_counts[1]*cell_counts[2];
	
	////scalar field property
	double* s=new double[vn];
	for(int i=0;i<vn;i++){
		s[i]=(double)i/(double)vn;
	}

	////vector field property
	double* v=new double[vn*3];
	for(int i=0;i<vn;i++){
		v[i*3]=.01;
		v[i*3+1]=v[i*3+2]=.0;
	}

	////IO test
	for(int frame=0;frame<100;frame++){
		px[3]+=.1;
		for(int i=0;i<vn;i++){
			v[i*3+1]+=.002;
		}

		std::string frame_dir=Frame_Dir(output_dir,frame);
		Create_Folder(output_dir,frame);

		{
			std::string file_name=frame_dir+"/grid";
			Write_Grid<double>(file_name,cell_counts,&dx,domain_min);
		}

		{
			std::string file_name=frame_dir+"/particles";
			Write_Particles<double>(file_name,px,pn);
		}	

		{
			std::string file_name=frame_dir+"/phi";
			Write_Scalar_Field<double>(file_name,s,cell_counts);
		}

		{
			std::string file_name=frame_dir+"/velocity";
			Write_Vector_Field<double>(file_name,v,cell_counts);
		}
	}
	
	delete [] v;
}

#endif