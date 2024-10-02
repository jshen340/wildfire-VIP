#pragma once
#ifdef USE_CPX

#include "Driver.h"
#include "Grid.h"
#include "Field.h"
#include "para.h"
#include "Poisson.h"
#include "CPXFunc.h"

class PoissonTestDriver : public Driver
{
public:
	int grid_size = 128;
	Grid<2> grid;

	Field<int, 2> cell_fixed;
	FaceField<Scalar, 2> face_vol;
	Field<Scalar, 2> cell_b;

	Poisson<2> poisson;

	Field<double, 2> field;

	virtual void Initialize()
	{
		grid.Initialize(Vector2i(grid_size, grid_size), 1.);

		poisson.Init(grid.cell_counts, 1.0);
		cell_fixed.Resize(Vector2i(grid_size, grid_size));
		face_vol.Resize(Vector2i(grid_size, grid_size));

		cell_b.Resize(Vector2i(grid_size, grid_size));

		init_boundary();

		field.Resize(Vector2i(grid_size, grid_size));
	}

	virtual void Write_Output_Files(int frame)
	{
		Driver::Write_Output_Files(frame);
		if (frame == 0)
		{
			std::string file_name = frame_dir + "/grid";
			grid.Write_To_File_3d(file_name);
			std::cout << "Write to file " << file_name << std::endl;
		}

		{
			std::string file_name = frame_dir + "/x";
			field.Write_To_File_3d(file_name);
		}
	}

	virtual void Run()
	{
		Timer timer;
		timer.Reset();

		timer.Elapse_And_Output_And_Reset("setup time");

		poisson.Solve();
		timer.Elapse_And_Output_And_Reset("solve time");

		for (int i = 0; i < grid_size; i++) for (int j = 0; j < grid_size; j++) field(Vector2i(i, j)) = (double)poisson.x(Vector2i(i, j));

		Write_Output_Files(0);
		timer.Elapse_And_Output_And_Reset("IO time");
	}

	void init_boundary()
	{
		cell_fixed.Fill(0);
		face_vol.Fill((Scalar)0);
		cell_b.Fill((Scalar)0);

		// all free surface
		if (test == 0)
		{
			face_vol.Fill((Scalar)1);

			cell_b.Fill((Scalar)1);
		}
		// wall on left
		if (test == 1)
		{
			face_vol.Fill((Scalar)1);
			for (int i = 0; i < grid_size; i++)
				face_vol(0, Vector2i(0, i)) = (Scalar)0;
			cell_b.Fill((Scalar)1);
		}
		// S1 domain
		if (test == 2)
		{
			face_vol.Fill((Scalar)1);

			int cx = grid_size / 2, cy = grid_size / 2;
			int cr = grid_size / 3;
			for (int i = 0; i < grid_size; i++) for (int j = 0; j < grid_size; j++)
			{
				if ((i - cx) * (i - cx) + (j - cy) * (j - cy) < cr * cr)
					cell_b(Vector2i(i, j)) = (Scalar)1;
				else
					cell_fixed(Vector2i(i, j)) = 1;
			}
		}
		// half wall inside
		if (test == 3)
		{
			face_vol.Fill((Scalar)1);
			for (int i = 0; i < grid_size / 2; i++)
				face_vol(0, Vector2i(grid_size / 3, i)) = (Scalar)0;
			cell_b.Fill((Scalar)1);
		}


		poisson.Init_Boundary(face_vol, cell_fixed, false);
		poisson.Update_b(cell_b);
	}
};

#endif