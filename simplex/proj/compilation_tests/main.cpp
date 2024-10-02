#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"


int main() {
	std::vector<glm::vec3> points;

	// generate points
	for (size_t i = 0; i < 3000; i++) {
		points.push_back(
			glm::vec3{ polyscope::randomUnit() - .5,
					  polyscope::randomUnit() - .5,
					  polyscope::randomUnit() - .5 });
	}

	std::vector<double> xC(points.size());
	for (size_t i = 0; i < points.size(); i++) {
		xC[i] = points[i].x;
	}
	

	// visualize!
	polyscope::init();
	polyscope::registerPointCloud("points", points);
	polyscope::getPointCloud("points")->addScalarQuantity("xC", xC);
	polyscope::show();
	return 0;
}