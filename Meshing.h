#pragma once
#include <iostream>
using namespace std;
#include "include/armadillo"
using namespace arma;

struct CNomalInfo
{
    
    CNomalInfo() = default;
    CNomalInfo(size_t etanum, size_t fainum)
    {
        eta_num = etanum;
        fai_num = fainum;
    }
    size_t eta_num{5};
    size_t fai_num{6};
};
struct CBlockInfo
{
    CBlockInfo() = default;
    CBlockInfo(double eta_delta, size_t eta_num, double fai_delta, size_t fai_num)
    {
        this->eta_delta = eta_delta;
        this->eta_num = eta_num;
        this->fai_delta = fai_delta;
        this->fai_num = fai_num;
    }
    double eta_delta{0.3}; // eta方向结构间隔
    size_t eta_num{4};     // eta方向每个block间网格数
    double fai_delta{0.3}; // fai方向结构间隔
    size_t fai_num{3};     // fai方向每个block间网格数
};

class Meshing
{
public:
    size_t eta_num{0};
    size_t fai_num{0};
    mat fai_mesh;
    mat eta_mesh;

private:
    void *m_mesh_res{nullptr};

public:
    Meshing() = default;
    ~Meshing();
    int calc_mesh(const CNomalInfo &info);
    int calc_mesh(const CBlockInfo &info);
private:
    void calc_mat();
};
