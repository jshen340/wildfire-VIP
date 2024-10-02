//////////////////////////////////////////////////////////////////////////
// Polyscope grid
// Copyright (c) (2018-), Bo Zhu, Zhecheng Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __PolyscopeGrid_h__
#define __PolyscopeGrid_h__

#include "MacGrid.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"
#include "PolyscopeViewerColor.h"

class PolyscopeGrid {
	Typedef_VectorDii(3);
public:
	std::string name;
	MacGrid<3> mac_grid;
	polyscope::CurveNetwork* vis_grid = nullptr;
	polyscope::PointCloud* vis_cell = nullptr;

	PolyscopeGrid() {}

	void Initialize(std::string grid_name, Grid<3>& _grid, bool enabled) {
		name = grid_name;
		mac_grid.Initialize(_grid);
		// draw grid
		Array<Vector3> nodes;
		Array<Vector2> edges;
		Grid<3>& grid = mac_grid.grid;
		int pair_counts = 0;
		//// x lines
		if (grid.node_counts[0] > 1) {
			for (int i = 0; i < grid.node_counts[1]; i++)for (int j = 0; j < grid.node_counts[2]; j++) {
				Vector3 start = grid.Node(Vector3i(0, i, j));
				Vector3 end = grid.Node(Vector3i(grid.node_counts[0] - 1, i, j));
				nodes.emplace_back(start);
				nodes.emplace_back(end);
				edges.emplace_back(2 * pair_counts, 2 * pair_counts + 1);
				pair_counts++;
			}
		}
		//// y lines
		if (grid.node_counts[1] > 1) {
			for (int i = 0; i < grid.node_counts[2]; i++)for (int j = 0; j < grid.node_counts[0]; j++) {
				Vector3 start = grid.Node(Vector3i(j, 0, i));
				Vector3 end = grid.Node(Vector3i(j, grid.node_counts[1] - 1, i));
				nodes.emplace_back(start);
				nodes.emplace_back(end);
				edges.emplace_back(2 * pair_counts, 2 * pair_counts + 1);
				pair_counts++;
			}
		}
		//// z lines
		if (grid.node_counts[2] > 1) {
			for (int i = 0; i < grid.node_counts[0]; i++)for (int j = 0; j < grid.node_counts[1]; j++) {
				Vector3 start = grid.Node(Vector3i(i, j, 0));
				Vector3 end = grid.Node(Vector3i(i, j, grid.node_counts[2] - 1)); nodes.emplace_back(start);
				nodes.emplace_back(end);
				edges.emplace_back(2 * pair_counts, 2 * pair_counts + 1);
				pair_counts++;
			}
		}
		vis_grid = polyscope::registerCurveNetwork(grid_name, nodes, edges);
		vis_grid->setColor(PolyscopeViewerColor::DimGrey);
		vis_grid->setMaterial("flat");
		vis_grid->setRadius(0.0005);
		vis_grid->setEnabled(enabled);
		// draw ghost cells
		Array<Vector3> cells;
		for (int i = 0; i < grid.cell_counts[0]; i++) {
			for (int j = 0; j < grid.cell_counts[1]; j++) {
				if (grid.cell_counts[2] == 0) {
					const Vector3i& cell = Vector3i(i, j, 0);
					const Vector3& pos = grid.Center(cell);
					cells.emplace_back(pos);
				} else {
					for (int k = 0; k < grid.cell_counts[2]; k++) {
						const Vector3i& cell = Vector3i(i, j, k);
						const Vector3& pos = grid.Center(cell);
						cells.emplace_back(pos);
					}
				}
			}
		}
		vis_cell = polyscope::registerPointCloud(grid_name + "_cells", cells);
		vis_cell->setPointColor(PolyscopeViewerColor::Black);
		vis_cell->setPointRenderMode(polyscope::PointRenderMode::Quad);
		vis_cell->setMaterial("flat");
		vis_cell->setTransparency(1.);
	}
};

#endif
