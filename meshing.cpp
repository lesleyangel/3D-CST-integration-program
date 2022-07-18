#include "Meshing.h"

extern "C"
{
    int print_and_2double(const char *);
    int addition(int, int);
    int sum_of_even(const int *, size_t);
    // test_struct handle_struct(test_struct);
    size_t count_char(const char *);
    // c_mesh_mat get_mesh_nomal(c_nomal_info);
    // Meshing::CMeshMat get_mesh_nomal(CNomalInfo);
    // Meshing::CMeshMat get_mesh_block(CBlockInfo);
    // void delete_mesh(Meshing::CMeshMat);

    void *get_mesh_nomal(CNomalInfo);
    void *get_mesh_block(CBlockInfo);
    void delete_mesh(void *);
    double get_eta(size_t, size_t, void *);
    double get_fai(size_t, size_t, void *);
    size_t get_eta_num(void *);
    size_t get_fai_num(void *);
    void print_license(void *);
}



int Meshing::calc_mesh(const CBlockInfo &info)
{
    m_mesh_res = get_mesh_block(info);
    calc_mat();
    return 0;
}

int Meshing::calc_mesh(const CNomalInfo &info)
{
    m_mesh_res = get_mesh_nomal(info);
    calc_mat();
    return 0;
}

void Meshing::calc_mat()
{
    print_license(m_mesh_res);
    fai_num = get_fai_num(m_mesh_res);
    eta_num = get_eta_num(m_mesh_res);
    fai_mesh.resize(fai_num, eta_num);
    eta_mesh.resize(fai_num, eta_num);
    for (size_t i = 0; i < fai_num; i++)
    {
        for (size_t j = 0; j < eta_num; j++)
        {
            fai_mesh(i, j) = get_fai(i, j, m_mesh_res);
            eta_mesh(i, j) = get_eta(i, j, m_mesh_res);
        }
    }
}

Meshing::~Meshing()
{
    if (m_mesh_res != nullptr)
    {
        delete_mesh(m_mesh_res);
        m_mesh_res = nullptr;
    }
}