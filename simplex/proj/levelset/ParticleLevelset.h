//#####################################################################
// Particle level set
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __ParticleLevelSet_h__
#define __ParticleLevelSet_h__
#include "LevelSet.h"
#include "Particles.h"
#include "RandomNumber.h"

template<int d> class ParticleLevelSet : public LevelSet<d>
{Typedef_VectorDii(d); using Base=LevelSet<d>;
public:
	using Base::grid; using Base::Phi; using Base::phi; using Base::Sign; using Base::Normal;
	Particles<d> particles;
	int particle_num_per_cell=Pow(4,d)*8;
	real band;
    real b_min;
    real b_max;
	real r_min;
	real r_max;
	// real target_phi_min;
	// real target_phi_max;
	RandomNumber random;
    
	//// data structure for reseeding
    Field<int,d> pcnt; ////the number of particles in each cell
    Field<int,d> grid_cell_to_array; //// map cell idx to array idx
    Array<int> array_to_grid_cell; //// map array idx to cell idx
    using PRI=std::pair<real,int>;
    Array<std::priority_queue<PRI,Array<PRI>,std::less<PRI> > > heaps; ////store the particles cell by cell in maxheap
	Array<Array<int>> arrays; ////store the particles cell by cell in array
    
	ParticleLevelSet(){}

	void Initialize(const Grid<d>& _grid)
	{
		Base::Initialize(_grid);
        pcnt.Resize(grid.cell_counts);
	    pcnt.Fill(0);
        grid_cell_to_array.Resize(grid.cell_counts);
        grid_cell_to_array.Fill(-1);
	}

	void Initialize_Particles()
	{
		band=grid.dx*(real)2;
		r_min=grid.dx*(real).1;
		r_max=grid.dx*(real).5;
        b_min=r_min;
        b_max=band;
		
		iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
			// VectorD pos=grid.Center(cell);	
			if(abs(phi(cell))<band){
				int e=particles.Add_Elements(particle_num_per_cell);
				for(int i=0;i<particle_num_per_cell;i++){
					VectorD pos=grid.dx*random.VectorValue<d>()+grid.Node(cell);
					particles.X(e+i)=pos;}}}
    
		
		for(int i=0;i<particles.Size();i++){
			// real target_phi=Sign(Phi(particles.X(i)))*random.Value()*(b_max-b_min)+b_min;
			
			bool valid=false;
			while(!valid){
				real target_phi=(random.Value()*(b_max-b_min)+b_min)*(i%2*2-1);
				real step=1;
				for(int iter=0;iter<15;iter++){
					particles.X(i)=particles.X(i)+(target_phi-Phi(particles.X(i)))*Normal(particles.X(i))*step;
					if(abs(Phi(particles.X(i)))<=b_max && abs(Phi(particles.X(i)))>=b_min) {valid=true; break;}
					else step*=0.5;
				}
				if(!valid)
					std::cout<<"Initialize_Particles invalid "<<i<<" "<<Phi(particles.X(i))<<" not in ["<<b_min<<","<<b_max<<"]"<<std::endl;
			}
			Particle_Sign(i)=Sign(Phi(particles.X(i)));
		}
		Update_Particle_Radii();
	}


	void Update_Particle_Radii()
	{
		for(int i=0;i<particles.Size();i++)
			Update_Single_Particle_Radii(i);
	}

	void Update_Single_Particle_Radii(int i){
		real p_phi=Phi(particles.X(i));real r_phi=abs(p_phi);
		// Particle_Sign(i)=Sign(p_phi);//// do not update sign here
		if(r_phi>r_max)Particle_Radius(i)=r_max;
		else if(r_phi<r_min)Particle_Radius(i)=r_min;
		else Particle_Radius(i)=r_phi;
		return;
		//// TOTRY: move particle ?
		if(r_phi<r_min*0.99){
			//// move towards normal
			particles.X(i)+=Normal(particles.X(i))*(r_min-r_phi)*Particle_Sign(i);
			p_phi=Phi(particles.X(i)); r_phi=abs(p_phi);
			Particle_Radius(i)=std::min(r_max, abs(r_phi));
			Particle_Sign(i)=p_phi<0?-1.0:1.0;
			// if(trap_in_loop>100)
			// 	std::cout<<"p_phi="<<p_phi<<",r_phi="<<r_phi<<",r_min="<<r_min<<"  "<<(r_phi<r_min*0.99)<<std::endl;
		}
		Particle_Sign(i)=p_phi<0?-1.0:1.0;
	}

	void Correction()
	{
		Field<real,d> phi_negative=phi;
		Field<real,d> phi_positive=phi;
		Grid<d> center_grid=Grid<d>(grid.cell_counts-VectorDi::Ones(),grid.dx,grid.domain_min+VectorD::Ones()*grid.dx*(real).5);
		for(int i=0;i<particles.Size();i++){
			real x_phi=Phi(particles.X(i));real p_phi=Particle_Radius(i);
			bool escaped=(Sign(x_phi)!=Particle_Sign(i))&&abs(x_phi)>p_phi;
			if(!escaped)continue; //// skip non-escaped particle

			VectorDi cell=center_grid.Cell_Coord(particles.X(i));
			for(int j=0;j<Grid<d>::Number_Of_Cell_Incident_Nodes();j++){
				VectorDi nb=Grid<d>::Cell_Incident_Node(cell,j);
				real node_phi=Particle_Phi(center_grid.Node(nb),i);
				if(Particle_Sign(i)>(real)0){phi_positive(nb)=std::max(phi_positive(nb),node_phi);}	////positive escaped particles
				else{phi_negative(nb)=std::min(phi_negative(nb),node_phi);}}}						////negative escaped particles

		iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
			phi(cell)=(abs(phi_positive(cell))<abs(phi_negative(cell)))?phi_positive(cell):phi_negative(cell);}
	}

	////
	void Reseed_Escaped_Particles(){ ////following "2003-Using the Particle Level Set Method and a Second Order Accurate Pressure Boundary Condition for Free Surface Flows"
		//// "In order to retain a smooth interface, escaped particles which are more than 1.5 times their radius removed from the appropriate side of the interface are deleted."
		std::set<int> remove_set;
		for(int i=0;i<particles.Size();i++){
			////remove invalid particles
			VectorDi cell=grid.Cell_Coord(particles.X(i));
			bool outside=!grid.Valid_Cell(cell);
			if(outside){remove_set.insert(i); continue;}
			//// remove far escaped particles
			real x_phi=Phi(particles.X(i));real p_phi=Particle_Radius(i);
			bool escaped=(Sign(x_phi)!=Particle_Sign(i))&&(abs(x_phi)>p_phi);
			if(!escaped)continue; //// skip non-escaped particle
			if(p_phi*1.5f<abs(x_phi)) remove_set.insert(i);}
		particles.Remove(remove_set);
		std::cout<<"Remove escaped particles:"<<remove_set.size()<<" "<<particles.Size()<<std::endl;
		remove_set.clear();
	}

	void Reseed(){ //// following "2002-A hybrid particle level set method for improved interface capturing"
        pcnt.Fill(0);
        grid_cell_to_array.Fill(-1);
        array_to_grid_cell.clear();
        heaps.clear();
        
		//// remove far p and obtain #p for cells
        std::set<int> remove_set;
		for(int i=0;i<particles.Size();i++){
			real x_phi=Phi(particles.X(i));real p_phi=Particle_Radius(i);
			bool escaped=(Sign(x_phi)!=Particle_Sign(i))&&abs(x_phi)>p_phi;
			if(escaped) continue; //// skip escaped particle

            if(x_phi<-band||x_phi>band) {remove_set.insert(i);} //// remove far particles
            else{VectorDi cell=grid.Cell_Coord(particles.X(i)); pcnt(cell)+=1;}} //// count close particles
		
		std::cout<<"Remove particles too far:"<<remove_set.size()<<" "<<particles.Size()<<std::endl;
        particles.Remove(remove_set);
        remove_set.clear();

        //// setup mapping between grid cell and array
        iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
            if(pcnt(cell)>0){
                grid_cell_to_array(cell)=array_to_grid_cell.size();
                array_to_grid_cell.push_back(grid.Cell_Index(cell));}}
        heaps.resize(array_to_grid_cell.size());

        //// build heap for each cell
        for(int i=0;i<particles.Size();i++){
            real x_phi=Phi(particles.X(i));real p_phi=Particle_Radius(i);
			bool escaped=(Sign(x_phi)!=Particle_Sign(i))&&abs(x_phi)>p_phi;
			if(escaped) continue;//// skip escaped particle

            VectorDi cell=grid.Cell_Coord(particles.X(i));
            if(grid_cell_to_array(cell)<0){std::cout<<"ERR: grid_cell_to_array(cell)<0 "<<cell.transpose()<<std::endl;
				std::cout<<i<<":"<<particles.X(i).transpose()<<std::endl;
				std::cout<<grid.cell_counts.transpose()<<std::endl;
			}
            heaps[grid_cell_to_array(cell)].push(PRI(abs(x_phi)-p_phi, i));}

		//// remove the non-escaped particles, if there are too many in a given cell 
        iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
            int cell_array_idx = grid_cell_to_array(cell);
            if(cell_array_idx<0) continue;

            auto &heap = heaps[grid_cell_to_array(cell)];
            //// remove particles in a cell
            while(heap.size()>1.6*particle_num_per_cell){
                // const real dist = heap.top().first;
                const int pidx = heap.top().second;
                heap.pop();
                remove_set.insert(pidx);}}
		
		std::cout<<"Remove particles in crowded cell:"<<remove_set.size()<<" "<<particles.Size()<<std::endl;
        particles.Remove(remove_set);
        remove_set.clear();

		//// insert new particles
		int insert_num;
        iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
            if(abs(phi(cell))>b_max) continue; ////skip far cell
            int o_cnt=pcnt(cell);
            if(o_cnt<=0) continue;
            if(o_cnt<particle_num_per_cell){//// insert particles
				insert_num+=particle_num_per_cell-o_cnt;
                int e=particles.Add_Elements(particle_num_per_cell-o_cnt);
                for(int i=e;i<particles.Size();i++){
					VectorD pos=grid.dx*random.VectorValue<d>()+grid.Node(cell);
					particles.X(i)=pos;
					
        			real target_phi=Phi(particles.X(i));
					if(abs(target_phi)<b_min) target_phi=Sign(target_phi)*b_min;
					if(abs(target_phi)>b_max) target_phi=Sign(target_phi)*b_max;
					
					while(!Attract_Particle(particles.X(i),target_phi,[&](VectorD p){
						real phi=Phi(p);
						return grid.Cell_Coord(p)==cell&&abs(phi)>b_min&&abs(phi)<b_max;
					})){
						particles.X(i)=grid.dx*random.VectorValue<d>()+grid.Node(cell);
					}

			        Particle_Sign(i)=Sign(Phi(particles.X(i)));
					Update_Single_Particle_Radii(i);}}}
		std::cout<<"Insert particles in empty cell:"<<insert_num<<" "<<particles.Size()<<std::endl;
	}

	bool Attract_Particle(VectorD& pos, real target_phi, std::function<bool(VectorD)> is_valid){
		bool valid=false;
		real step=1;
		VectorD p=pos;
		for(int iter=0;iter<15;iter++){
			p=p+(target_phi-Phi(p))*Normal(p)*step;
			if(is_valid(pos)){valid=true; break;}
			else step*=0.5;
		}
		// if(!valid)
		// 	std::cout<<"Attract_Particle invalid "<<pos.transpose()<<"->"<<p.transpose()<<","<<Phi(p)<<" invalid"<<std::endl;
		pos=p;
		return valid;
	}

	void Reseed_Random(){
 		pcnt.Fill(0);
        grid_cell_to_array.Fill(-1);
        array_to_grid_cell.clear();
        arrays.clear();
        
        std::set<int> remove_set;
		for(int i=0;i<particles.Size();i++){
			real x_phi=Phi(particles.X(i));real p_phi=Particle_Radius(i);
			bool escaped=(Sign(x_phi)!=Particle_Sign(i))&&abs(x_phi)>p_phi;
			if(escaped) continue; //// skip escaped particle

            if(x_phi<-b_max||x_phi>b_max) {remove_set.insert(i);} //// remove far particles
            else{VectorDi cell=grid.Cell_Coord(particles.X(i)); pcnt(cell)+=1;}} //// count close particles
        particles.Remove(remove_set);
        remove_set.clear();

        //// setup mapping between grid cell and array
        iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
            if(pcnt(cell)>0){
                grid_cell_to_array(cell)=array_to_grid_cell.size();
                array_to_grid_cell.push_back(grid.Cell_Index(cell));}}
        arrays.resize(array_to_grid_cell.size());

        //// build array for each cell
        for(int i=0;i<particles.Size();i++){
            real x_phi=Phi(particles.X(i));real p_phi=Particle_Radius(i);
			bool escaped=(Sign(x_phi)!=Particle_Sign(i))&&abs(x_phi)>p_phi;
			if(escaped) continue;//// skip escaped particle

            VectorDi cell=grid.Cell_Coord(particles.X(i));
            if(grid_cell_to_array(cell)<0){std::cout<<"ERR: grid_cell_to_array(cell)<0"<<cell.transpose()<<std::endl;
				std::cout<<particles.X(i).transpose()<<std::endl;
				std::cout<<grid.cell_counts.transpose()<<std::endl;}
			arrays[grid_cell_to_array(cell)].push_back(i);}

		//// remove the non-escaped particles, if there are too many in a given cell 
        iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
            int cell_array_idx = grid_cell_to_array(cell);
            if(cell_array_idx<0) continue;

            auto &array = arrays[grid_cell_to_array(cell)];
            //// remove particles in a cell
			if(array.size()>1.6*particle_num_per_cell){
				int remove_size=array.size()-1.6*particle_num_per_cell;
				while(remove_set.size()<remove_size){
					remove_set.insert(array[int(random.Value()*(array.size()-1))]);}
				particles.Remove(remove_set);
        		remove_set.clear();}}
                
        iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
            if(abs(phi(cell))>band) continue; ////skip far cell
            int o_cnt=pcnt(cell);
            if(o_cnt<=0) continue;
            if(o_cnt<particle_num_per_cell){//// insert particles
                int e=particles.Add_Elements(particle_num_per_cell-o_cnt);
                for(int i=0;i<particle_num_per_cell-o_cnt;i++){
					VectorD pos=grid.dx*random.VectorValue<d>()+grid.Node(cell);
					particles.X(e+i)=pos;
        			// real target_phi=(random.Value()*(b_max-b_min)-b_min)*(i%2*2-1);
					real target_phi=std::min(abs(Phi(particles.X(i))), b_max);
					target_phi=std::max(target_phi, b_min)*Sign(Phi(particles.X(i)));
					
					real step=1;
					for(int iter=0;iter<5;iter++){
						particles.X(e+i)=pos+(target_phi-Phi(pos))*Normal(pos)*step;
						if(abs(Phi(particles.X(i)))<band) break;
						else step*=0.5;}
			        Particle_Sign(e+i)=Sign(Phi(particles.X(e+i)));
					Update_Single_Particle_Radii(e+i);}}}
	}

	

public:
	real& Particle_Radius(const int idx){return particles.C(idx);}
	const real& Particle_Radius(const int idx) const {return particles.C(idx);}
	real& Particle_Sign(const int idx){return particles.M(idx);}
	const real& Particle_Sign(const int idx) const {return particles.M(idx);}
	real Particle_Phi(const VectorD& pos,const int idx) const {return particles.M(idx)*(Particle_Radius(idx)-(pos-particles.X(idx)).norm());}
};

#endif