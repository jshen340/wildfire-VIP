#include "randomTest.h"
#include "form01mapping.h"
#include "cuda_runtime_api.h"
#include "gpuUtils.h"
#include "cublas_v2.h"
#include <memory>

void printFace(const grid2D& grid, Scalar *v)
{
	printf("--------------\n");
	for (int j = 0; j < grid.Ny; j++)
	{
		for (int i = 0; i <= grid.Nx; i++)
			printf("%f ", v[grid.face_ind(i, j, 0)]);
		printf("\n");
	}

	printf("\t----\t\n");

	for (int j = 0; j <= grid.Ny; j++)
	{
		for (int i = 0; i < grid.Nx; i++)
			printf("%f ", v[grid.face_size(0) + grid.face_ind(i, j, 1)]);
		printf("\n");
	}
	printf("--------------\n");
}

void printNode(const grid2D& grid, Scalar *v)
{
	printf("--------------\n");
	for (int j = 0; j <= grid.Ny; j++)
	{
		for (int i = 0; i <= grid.Nx; i++)
			printf("%f ", v[grid.node_ind(i, j)]);
		printf("\n");
	}
	printf("--------------\n");
}

void curlTest()
{
	int N = 16;
	grid2D grid;
	grid.init(N, N);

	Scalar *h_grad, *d_grad;
	h_grad = (Scalar*)malloc(sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));
	cudaMalloc(&d_grad, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));

	memset(h_grad, 0, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));

	for (int i = 0; i < N-1; i++)
		for (int j = 0; j < N; j++)
			h_grad[grid.face_ind(i + 1, j, 0)] = (i + 1)*j - i * j;

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N-1; j++)
			h_grad[grid.face_size(0) + grid.face_ind(i, j + 1, 1)] = i * (j + 1) - i * j;

	printFace(grid, h_grad);

	cudaMemcpy(d_grad, h_grad, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)), cudaMemcpyHostToDevice);

	Scalar *h_curl, *d_curl;
	h_curl = (Scalar*)malloc(sizeof(Scalar)*grid.node_size());
	cudaMalloc(&d_curl, sizeof(Scalar)*grid.node_size());

	grid2DOperator::Cod1Mapping(grid, d_curl, d_grad, d_grad + grid.face_size(0));

	cudaMemcpy(h_curl, d_curl, sizeof(Scalar)*grid.node_size(), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	printNode(grid, h_curl);
}

void gradTest()
{
	int N = 16;
	grid2D grid;
	grid.init(N, N);

	Scalar *h_curl, *d_curl;
	h_curl = (Scalar*)malloc(sizeof(Scalar)*grid.node_size());
	cudaMalloc(&d_curl, sizeof(Scalar)*grid.node_size());


	memset(h_curl, 0, sizeof(Scalar)*grid.node_size());

	for (int i = 0; i <= N; i++)
		for (int j = 0; j <= N; j++)
			h_curl[grid.node_ind(i, j)] = i * j;

	printNode(grid, h_curl);

	cudaMemcpy(d_curl, h_curl, sizeof(Scalar)*grid.node_size(), cudaMemcpyHostToDevice);

	Scalar *h_grad, *d_grad;
	h_grad = (Scalar*)malloc(sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));
	cudaMalloc(&d_grad, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));

	grid2DOperator::D0Mapping(grid, d_grad, d_grad + grid.face_size(0), d_curl);

	cudaMemcpy(h_grad, d_grad, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	printFace(grid, h_grad);
}

void vecLapTest()
{
	int N = 8;
	grid2D grid;
	grid.init(N, N);

	Scalar *h_v, *d_v;
	h_v = (Scalar*)malloc(sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));
	cudaMalloc(&d_v, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));

	memset(h_v, 0, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));

	for (int i = 0; i <= N; i++)
		for (int j = 0; j < N; j++)
			h_v[grid.face_ind(i, j, 0)] = i * i;

	for (int i = 0; i < N; i++)
		for (int j = 0; j <= N; j++)
			h_v[grid.face_size(0) + grid.face_ind(i, j, 1)] = j * j;

	printFace(grid, h_v);

	cudaMemcpy(d_v, h_v, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)), cudaMemcpyHostToDevice);

	Scalar *h_temp, *d_temp;
	h_temp = (Scalar*)malloc(sizeof(Scalar)*(grid.cell_size() + grid.face_size(0) + grid.face_size(1) + grid.node_size()));
	cudaMalloc(&d_temp, sizeof(Scalar)*(grid.cell_size() + grid.face_size(0) + grid.face_size(1) + grid.node_size()));

	Scalar *h_lap, *d_lap;
	h_lap = (Scalar*)malloc(sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));
	cudaMalloc(&d_lap, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));

	cudaMemset(d_lap, 0, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));

	auto neg = [=]__device__(Scalar& tv) { tv = -tv; };

	int cell_off = 0;
	int face_x_off = cell_off + grid.cell_size();
	int face_y_off = face_x_off + grid.face_size(0);
	int node_off = face_y_off + grid.face_size(1);
	Scalar mu = (Scalar)1;
	cublasHandle_t cublasHandle;
	cublasCreate(&cublasHandle);

	// grad div v
	cudaMemcpy(d_temp + face_x_off, d_v, sizeof(Scalar)*grid.face_size(0), cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_temp + face_y_off, d_v + grid.face_size(0), sizeof(Scalar)*grid.face_size(1), cudaMemcpyDeviceToDevice);

	cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));

	grid2DOperator::D1Mapping(grid, d_temp + cell_off, d_temp + face_x_off, d_temp + face_y_off);

	grid2DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + cell_off);

	Axpy(cublasHandle, grid.face_size(0), &mu, d_temp + face_x_off, 1, d_lap, 1);
	Axpy(cublasHandle, grid.face_size(1), &mu, d_temp + face_y_off, 1, d_lap + grid.face_size(0), 1);

	// - curl curl v
	grid2DOperator::Cod1Mapping(grid, d_temp + cell_off, d_v, d_v +grid.face_size(0));

	grid2DOperator::D0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + cell_off);

	//cwise_mapping_wrapper(d_temp + face_x_off, neg, grid.face_size(0));
	cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));

	Axpy(cublasHandle, grid.face_size(0), &mu, d_temp + face_x_off, 1, d_lap, 1);
	Axpy(cublasHandle, grid.face_size(1), &mu, d_temp + face_y_off, 1, d_lap + grid.face_size(0), 1);

	cudaMemcpy(h_lap, d_lap, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	printFace(grid, h_lap);
}