#ifndef __IO_Adapter_h__
#define __IO_Adapter_h__
#include <fstream>
#include <iostream>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <string>
#include <vector>
#ifdef WIN32
#include <windows.h>
#elif defined(__linux__)
#include <sys/stat.h>
#endif

#ifdef WIN32
#define NO_MINMAX
#undef min
#undef max
#endif

namespace IOAdapter
{
template<class T_VAL> void Write_Binary(std::ostream& output,const T_VAL& data)
{output.write(reinterpret_cast<const char*>(&data),sizeof(T_VAL));}

template<class T_VAL> void Write_Binary_Array(std::ostream& output,const T_VAL* array,const int n)
{if(n>0)output.write(reinterpret_cast<const char*>(array),n*sizeof(T_VAL));}

template<class T_VAL> bool Write_Binary_To_File(const std::string& file_name,const T_VAL& data)
{std::ofstream output(file_name,std::ios::binary);if(!output)return false;Write_Binary(output,data);return true;}

template<class T_VAL> bool Write_Binary_Array_To_File(const std::string& file_name,T_VAL* array,const int n)
{std::ofstream output(file_name,std::ios::binary);if(!output)return false;Write_Binary_Array(output,array,n);return true;}

template<class T_VAL> void Write_Text(std::ostream& output,const T_VAL& data)
{output<<data;}

template<class T_VAL> bool Write_Text_To_File(const std::string& file_name,const T_VAL& data)
{std::ofstream output(file_name);if(!output)return false;Write_Text(output,data);return true;}

#ifdef WIN32
inline bool Directory_Exists(const char* dirname)
{DWORD attr=GetFileAttributes(dirname);return((attr!=-1)&&(attr&FILE_ATTRIBUTE_DIRECTORY));}

inline bool Create_Directory(const std::string& dirname)
{
    if(!Directory_Exists(dirname.c_str())){size_t pos=0;
        do{pos=dirname.find_first_of("\\/",pos+1);
        if(!Directory_Exists(dirname.substr(0,pos).c_str())){
            if(CreateDirectory(dirname.substr(0,pos).c_str(),NULL)==0 && ERROR_ALREADY_EXISTS!=GetLastError()){
                std::cerr<<"Error: [File] Create directory "<<dirname<<" failed!"<<std::endl;return false;}}}while(pos!=std::string::npos);}
    return true;
}
#elif defined(__linux__)
inline bool Directory_Exists(const char* dirname)
{struct stat s;return stat(dirname,&s)==0;}

inline bool Create_Directory(const std::string& dirname)
{if(!Directory_Exists(dirname.c_str())){size_t pos=0;
        do{pos=dirname.find_first_of("\\/",pos+1);
            if(!Directory_Exists(dirname.substr(0,pos).c_str())){
                if(mkdir(dirname.substr(0,pos).c_str(),0755)!=0 && errno!=EEXIST){
                    std::cerr<<"Error: [File] Create directory "<<dirname<<"failed!"<<std::endl;return false;}}}
		while(pos!=std::string::npos);}
return true;}
#endif

//////////////////////////////////////////////////////////////////////////
////IO functions
//////////////////////////////////////////////////////////////////////////

std::string Frame_Dir(const std::string& output_dir,const int frame)
{
	return output_dir+"/"+std::to_string(frame);
}

void Create_Folder(const std::string& output_dir,const int frame)
{
	if(!Directory_Exists(output_dir.c_str()))Create_Directory(output_dir);

	std::string frame_dir=Frame_Dir(output_dir,frame);
	if(!Directory_Exists(frame_dir.c_str()))Create_Directory(frame_dir);
		
	{std::string file_name=output_dir+"/0/last_frame.txt";
	Write_Text_To_File(file_name,std::to_string(frame));}
}

template<class T> void Write_Grid(const std::string& file_name,const int* cell_counts,const T* dx,const T* domain_min)
{
    std::ofstream output(file_name,std::ios::binary);
	if(!output){std::cerr<<"Cannot open file "<<file_name<<std::endl;return;}

	Write_Binary_Array(output,cell_counts,3);
    Write_Binary(output,*dx);
    Write_Binary_Array(output,domain_min,3);
}

template<class T> void Write_Points(const std::string& file_name,const T* x,const int n)
{
	std::ofstream output(file_name,std::ios::binary);
	if(!output){std::cerr<<"Cannot open file "<<file_name<<std::endl;return;}

	Write_Binary(output,n);
	if(n>0){Write_Binary_Array(output,x,n*3);}
}

template<class T> void Write_Particles(const std::string& file_name,const T* x,const int n,const T* v=nullptr,const T* f=nullptr,const T* m=nullptr)
{
	std::ofstream output(file_name,std::ios::binary);
	if(!output){std::cerr<<"Cannot open file "<<file_name<<std::endl;return;}

	Write_Binary(output,n);
	if(n==0)return;

	std::vector<T> placeholder_T(n*3,(T)0);
	std::vector<int> placeholder_i(n,0);

	if(n>0){Write_Binary_Array(output,x,n*3);}

	if(v)Write_Binary_Array(output,v,n*3);
	else Write_Binary_Array(output,&placeholder_T[0],n*3);

	if(f)Write_Binary_Array(output,f,n*3);
	else Write_Binary_Array(output,&placeholder_T[0],n*3);

	if(m)Write_Binary_Array(output,m,n);
	else Write_Binary_Array(output,&placeholder_T[0],n);

	Write_Binary_Array(output,&placeholder_T[0],n);

	Write_Binary_Array(output,&placeholder_i[0],n);
}

template<class T> void Write_Scalar_Field(const std::string& file_name,const T* s,const int* counts)
{
	std::ofstream output(file_name,std::ios::binary);
	if(!output){std::cerr<<"Cannot open file "<<file_name<<std::endl;return;}

	Write_Binary_Array(output,counts,3);

	int n=counts[0]*counts[1]*counts[2];
	Write_Binary_Array(output,s,n);
}

template<class T> void Write_Vector_Field(const std::string& file_name,const T* v,const int* counts)
{
	std::ofstream output(file_name,std::ios::binary);
	if(!output){std::cerr<<"Cannot open file "<<file_name<<std::endl;return;}

	Write_Binary_Array(output,counts,3);

	int n=counts[0]*counts[1]*counts[2]*3;
	Write_Binary_Array(output,v,n);
}

};

#endif