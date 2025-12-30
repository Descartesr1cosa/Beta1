#pragma once
#include "3_field/Halo_Type.h"
#include "3_field/2_MPCNS_Field.h"
#include "0_basic/MPI_WRAPPER.h"

class Halo
{
public:
    Halo(Field *field, TOPO::Topology *topo)
    {
        fld_ = field;
        topo_ = topo;
        build_pattern();
    };

    void data_trans(std::string &field_name)
    {
        exchange_inner(field_name);
        exchange_parallel(field_name);
    }

    void data_trans_2DCorner(std::string &field_name)
    {
        exchange_inner_edge(field_name);
        exchange_parallel_edge(field_name);
    }

    void data_trans_3DCorner(std::string &field_name)
    {
        exchange_inner_vertex(field_name);
        exchange_parallel_vertex(field_name);
    }

    // face
    void exchange_inner(std::string field_name);
    void exchange_parallel(std::string field_name);

    // edge
    void exchange_inner_edge(std::string field_name);
    void exchange_parallel_edge(std::string field_name);

    // vertex
    void exchange_inner_vertex(std::string field_name);
    void exchange_parallel_vertex(std::string field_name);

private:
    Field *fld_;
    TOPO::Topology *topo_;

    // key = (StaggerLocation, nghost)，value = HaloPattern
    using PatternKey = std::pair<StaggerLocation, int>;
    std::map<PatternKey, HaloPattern> inner_patterns_;
    std::map<PatternKey, HaloPattern> parallel_patterns_;

    // For 2D Corner
    std::map<PatternKey, HaloPattern> inner_edge_patterns_;
    std::map<PatternKey, HaloPattern> parallel_edge_patterns_send;
    std::map<PatternKey, HaloPattern> parallel_edge_patterns_recv;

    // For 3D Corner
    std::map<PatternKey, HaloPattern> inner_vertex_patterns_;
    std::map<PatternKey, HaloPattern> parallel_vertex_patterns_send;
    std::map<PatternKey, HaloPattern> parallel_vertex_patterns_recv;

    // 复用的 send / recv 缓冲区（MPI 并行用）
    std::vector<std::vector<double>> send_buf;
    std::vector<std::vector<double>> recv_buf;
    std::vector<MPI_Request> req_send, req_recv;
    std::vector<MPI_Status> stat_send, stat_recv;
    std::vector<int32_t> length, length_corner_recv;

    //=========================================================================

    void build_pattern();

    void build_inner_pattern(StaggerLocation loc, int nghost);
    void build_parallel_pattern(StaggerLocation loc, int nghost);
    void build_inner_edge_pattern(StaggerLocation loc, int nghost);
    void build_parallel_edge_pattern(StaggerLocation loc, int nghost);
    void build_inner_vertex_pattern(StaggerLocation loc, int nghost);
    void build_parallel_vertex_pattern(StaggerLocation loc, int nghost);

    void mpi_exchange_edge_meta(
        const std::map<int, std::vector<EdgeMeta>> &meta_to_send,
        std::vector<EdgeMeta> &recv_metas);
    void mpi_exchange_vertex_meta(
        const std::map<int, std::vector<VertexMeta>> &meta_to_send,
        std::vector<VertexMeta> &recv_metas);

    void switch_range_ghost(Box3 &box, Direction &dir, int nghost);
    void switch_range_inner(Box3 &box, Direction &dir, int nghost);

    Box3 make_cell_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_cell_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facex_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facex_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facey_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facey_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facez_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_facez_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgex_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgex_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgey_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgey_inner_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgez_ghost_box(Box3 &face, Direction &dir, int nghost);
    Box3 make_edgez_inner_box(Box3 &face, Direction &dir, int nghost);

    // 在 Cell 空间里，为二维角区（edge strip：两个方向 ghost）构造 box
    // edge_node 是该 edge 在 node 空间里的交集区域（一般是一条线段 [lo,hi)）
    // 接受可以直接两个方向均为ghost，代表角区
    // 发送则一定是一个为inner一个为ghost，约定，dir1为inner
    Box3 make_cell_edge_ghost_box(Box3 &edge_node, Direction &dir1, Direction &dir2, int nghost);
    Box3 make_cell_edge_innerghost_box(Box3 &edge_node, Direction &dir1, Direction &dir2, int nghost);

    // 在 Cell 空间里，为三维角区（vertex：三个方向 ghost）构造 box
    // vertex_node 是顶点附近的 node 区域（通常 [i,i+1)×[j,j+1)×[k,k+1)）
    Box3 make_cell_vertex_ghost_box(Box3 &vertex_node, Direction &dir1, Direction &dir2, Direction &dir3, int nghost);
    Box3 make_cell_vertex_innerghost_box(Box3 &vertex_node, Direction &dir1, Direction &dir2, Direction &dir3, int nghost);
};
