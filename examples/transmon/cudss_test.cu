// cuDSS sparse direct solver test
// Reads Matrix Market file and binary RHS, solves with cuDSS
//
// Compile:
//   nvcc cudss_test.cu -I${CUDSS_DIR}/include -L${CUDSS_DIR}/lib -lcudss -lcublas -o cudss_test

#include <cudss.h>
#include <cuda_runtime.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>

#define CHECK_CUDA(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl; \
        exit(1); \
    } \
} while(0)

#define CHECK_CUDSS(call) do { \
    cudssStatus_t status = call; \
    if (status != CUDSS_STATUS_SUCCESS) { \
        std::cerr << "cuDSS error: " << status << std::endl; \
        exit(1); \
    } \
} while(0)

void readMatrixMarket(const std::string &filename,
                      int64_t &nrows, int64_t &ncols, int64_t &nnz,
                      std::vector<int> &rowPtr,
                      std::vector<int> &colIdx,
                      std::vector<double> &values)
{
    std::ifstream file(filename);
    if (!file) { std::cerr << "Cannot open: " << filename << std::endl; exit(1); }

    std::string line;
    while (std::getline(file, line) && line[0] == '%');

    std::istringstream(line) >> nrows >> ncols >> nnz;

    std::vector<int> cooRow(nnz), cooCol(nnz);
    std::vector<double> cooVal(nnz);
    for (int64_t i = 0; i < nnz; i++) {
        int r, c; double v;
        file >> r >> c >> v;
        cooRow[i] = r - 1; cooCol[i] = c - 1; cooVal[i] = v;
    }

    rowPtr.resize(nrows + 1, 0);
    for (int64_t i = 0; i < nnz; i++) rowPtr[cooRow[i] + 1]++;
    for (int64_t i = 0; i < nrows; i++) rowPtr[i + 1] += rowPtr[i];

    colIdx.resize(nnz); values.resize(nnz);
    std::vector<int> cnt(nrows, 0);
    for (int64_t i = 0; i < nnz; i++) {
        int row = cooRow[i], dest = rowPtr[row] + cnt[row]++;
        colIdx[dest] = cooCol[i]; values[dest] = cooVal[i];
    }
}

// Read binary vector file (int size, then doubles)
void readBinaryVector(const std::string &filename, std::vector<double> &vec)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file) { std::cerr << "Cannot open: " << filename << std::endl; exit(1); }
    int size;
    file.read((char*)&size, sizeof(int));
    vec.resize(size);
    file.read((char*)vec.data(), size * sizeof(double));
}

double computeResidual(int64_t n, const std::vector<int> &rowPtr,
                       const std::vector<int> &colIdx, const std::vector<double> &val,
                       const std::vector<double> &x, const std::vector<double> &b)
{
    double res = 0, bnorm = 0;
    for (int64_t i = 0; i < n; i++) {
        double Axi = 0;
        for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) Axi += val[j] * x[colIdx[j]];
        res += (Axi - b[i]) * (Axi - b[i]);
        bnorm += b[i] * b[i];
    }
    return std::sqrt(res / bnorm);
}

int main(int argc, char *argv[])
{
    std::string matFile = argc > 1 ? argv[1] : "matrix.mtx";
    std::string rhsFile = argc > 2 ? argv[2] : "rhs.bin";
    std::string solFile = argc > 3 ? argv[3] : "sol.bin";
    int numGpus = argc > 4 ? std::atoi(argv[4]) : 1;
    int ndLevels = argc > 5 ? std::atoi(argv[5]) : 0;  // 0 = default

    int availGpus;
    CHECK_CUDA(cudaGetDeviceCount(&availGpus));
    numGpus = std::min(numGpus, availGpus);
    std::cout << "Using " << numGpus << " GPU(s)";
    if (ndLevels > 0) std::cout << ", ND levels: " << ndLevels;
    std::cout << "\n";

    int64_t n, nc, nnz;
    std::vector<int> rowPtr, colIdx;
    std::vector<double> val, rhs;
    readMatrixMarket(matFile, n, nc, nnz, rowPtr, colIdx, val);
    std::cout << "Matrix: " << n << " x " << nc << ", nnz = " << nnz << std::endl;
    
    readBinaryVector(rhsFile, rhs);

    int *d_row, *d_col; double *d_val, *d_rhs, *d_sol;
    CHECK_CUDA(cudaMalloc(&d_row, (n+1)*sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_col, nnz*sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_val, nnz*sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_rhs, n*sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_sol, n*sizeof(double)));

    CHECK_CUDA(cudaMemcpy(d_row, rowPtr.data(), (n+1)*sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_col, colIdx.data(), nnz*sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_val, val.data(), nnz*sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_rhs, rhs.data(), n*sizeof(double), cudaMemcpyHostToDevice));

    cudssHandle_t handle; cudssConfig_t config; cudssData_t data;
    cudssMatrix_t A, B, X;

    if (numGpus == 1) {
        CHECK_CUDSS(cudssCreate(&handle));
    } else {
        std::vector<int> deviceIds(numGpus);
        for (int i = 0; i < numGpus; i++) deviceIds[i] = i;
        CHECK_CUDSS(cudssCreateMg(&handle, numGpus, deviceIds.data()));
    }

    // Enable MT mode with OpenMP threading layer (set CUDSS_THREADING_LIB env var)
    const char* threadLib = std::getenv("CUDSS_THREADING_LIB");
    if (threadLib) {
        CHECK_CUDSS(cudssSetThreadingLayer(handle, threadLib));
        std::cout << "Using threading layer: " << threadLib << "\n";
    }
    CHECK_CUDSS(cudssConfigCreate(&config));
    
    // Set ND levels if specified
    if (ndLevels > 0) {
        cudssConfigSet(config, CUDSS_CONFIG_ND_NLEVELS, &ndLevels, sizeof(ndLevels));
    }
    
    // Try hybrid mode for multi-GPU (can be faster for analysis)
    if (numGpus > 1) {
        int hybrid = 1;
        cudssConfigSet(config, CUDSS_CONFIG_HYBRID_MODE, &hybrid, sizeof(hybrid));
    }
    
    CHECK_CUDSS(cudssDataCreate(handle, &data));
    CHECK_CUDSS(cudssMatrixCreateCsr(&A, n, nc, nnz, d_row, nullptr, d_col, d_val,
        CUDA_R_32I, CUDA_R_64F, CUDSS_MTYPE_GENERAL, CUDSS_MVIEW_FULL, CUDSS_BASE_ZERO));
    CHECK_CUDSS(cudssMatrixCreateDn(&B, n, 1, n, d_rhs, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR));
    CHECK_CUDSS(cudssMatrixCreateDn(&X, n, 1, n, d_sol, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR));

    auto t0 = std::chrono::high_resolution_clock::now();
    CHECK_CUDSS(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, config, data, A, X, B));
    CHECK_CUDA(cudaDeviceSynchronize());
    auto t1 = std::chrono::high_resolution_clock::now();
    CHECK_CUDSS(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, config, data, A, X, B));
    CHECK_CUDA(cudaDeviceSynchronize());
    auto t2 = std::chrono::high_resolution_clock::now();
    CHECK_CUDSS(cudssExecute(handle, CUDSS_PHASE_SOLVE, config, data, A, X, B));
    CHECK_CUDA(cudaDeviceSynchronize());
    auto t3 = std::chrono::high_resolution_clock::now();

    std::cout << "Analysis:      " << std::chrono::duration<double>(t1-t0).count() << " s\n";
    std::cout << "Factorization: " << std::chrono::duration<double>(t2-t1).count() << " s\n";
    std::cout << "Solve:         " << std::chrono::duration<double>(t3-t2).count() << " s\n";
    std::cout << "Total:         " << std::chrono::duration<double>(t3-t0).count() << " s\n";

    std::vector<double> sol(n);
    CHECK_CUDA(cudaMemcpy(sol.data(), d_sol, n*sizeof(double), cudaMemcpyDeviceToHost));

    std::cout << "Residual:      " << computeResidual(n, rowPtr, colIdx, val, sol, rhs) << std::endl;

    // Compare with Palace reference solution
    std::vector<double> refSol;
    readBinaryVector(solFile, refSol);
    
    if (refSol.size() == sol.size()) {
        double diff = 0, refNorm = 0;
        for (size_t i = 0; i < sol.size(); i++) {
            diff += (sol[i] - refSol[i]) * (sol[i] - refSol[i]);
            refNorm += refSol[i] * refSol[i];
        }
        std::cout << "||cuDSS - Palace|| / ||Palace|| = " << std::sqrt(diff/refNorm) << std::endl;
    }

    CHECK_CUDSS(cudssMatrixDestroy(A));
    CHECK_CUDSS(cudssMatrixDestroy(B));
    CHECK_CUDSS(cudssMatrixDestroy(X));
    CHECK_CUDSS(cudssDataDestroy(handle, data));
    CHECK_CUDSS(cudssConfigDestroy(config));
    CHECK_CUDSS(cudssDestroy(handle));
    cudaFree(d_row); cudaFree(d_col); cudaFree(d_val); cudaFree(d_rhs); cudaFree(d_sol);

    return 0;
}
