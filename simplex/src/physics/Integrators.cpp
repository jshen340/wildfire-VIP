#include "Integrators.h"

template class IntegratorX<2>;
template class IntegratorX<3>;

template<int d>
void IntegratorXImplicit<d>::Integrate(Array<VectorD>& X, Array<VectorD>& V, const real dt)
{
	int n = X.size();
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		X[i] += V[i] * dt;
		real dis; VectorD normal;
		if (boundary.Nearest_Boundary(X[i], dis, normal)) {
			const VectorD& v = V[i];
			real vn_rel = normal.dot(v);
			if (dis <= 0) {
				X[i] -= normal * dis;
				VectorD vn = vn_rel * normal;
				VectorD vt = v - vn;
				V[i] = vt * t_coeff + vn * n_coeff;
			}
		}
		//if (vn_rel * dt + dis <= check_width);
	}
}

template class IntegratorXImplicit<2>;
template class IntegratorXImplicit<3>;