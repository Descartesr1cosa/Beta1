#include "3_field/3_MPCNS_Halo.h"

void Halo::build_inner_vertex_pattern(StaggerLocation loc, int nghost)
{

    auto block_node_size = [](const Block &blk) -> Int3
    {
        return {blk.mx, blk.my, blk.mz};
    };

    auto detect_direction = [](Box3 &face, Int3 &blk_mxyz, int &dir1, int &dir2, int &dir3)
    {
        dir1 = -999;
        dir2 = -999;
        dir3 = -999;
        // 下面的判断逻辑基于假设：
        // 接口在某一方向上厚度为 1 层节点，而且刚好贴在 0 或 Ni_node 上。

        // XMinus: i=0 这一层节点
        if (face.lo.i == 0 && face.hi.i == 1)
        {
            int direction = -1;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
                dir2 = direction;
            else if (dir3 == -999)
            {
                dir3 = direction;
                return;
            }
        }

        // XPlus: i=Ni_node 这一层节点
        if (face.lo.i == blk_mxyz.i && face.hi.i == blk_mxyz.i + 1)
        {
            int direction = 1;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
                dir2 = direction;
            else if (dir3 == -999)
            {
                dir3 = direction;
                return;
            }
        }

        // YMinus: j=0
        if (face.lo.j == 0 && face.hi.j == 1)
        {
            int direction = -2;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
                dir2 = direction;
            else if (dir3 == -999)
            {
                dir3 = direction;
                return;
            }
        }

        // YPlus: j=Nj_node
        if (face.lo.j == blk_mxyz.j && face.hi.j == blk_mxyz.j + 1)
        {
            int direction = 2;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
                dir2 = direction;
            else if (dir3 == -999)
            {
                dir3 = direction;
                return;
            }
        }

        // ZMinus: k=0
        if (face.lo.k == 0 && face.hi.k == 1)
        {
            int direction = -3;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
                dir2 = direction;
            else if (dir3 == -999)
            {
                dir3 = direction;
                return;
            }
        }

        // ZPlus: k=Nk_node
        if (face.lo.k == blk_mxyz.k && face.hi.k == blk_mxyz.k + 1)
        {
            int direction = 3;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
                dir2 = direction;
            else if (dir3 == -999)
            {
                dir3 = direction;
                return;
            }
        }

        std::cout << "block:" << blk_mxyz.i << " " << blk_mxyz.j << " " << blk_mxyz.k << "\n"
                  << face.lo.i << "\t" << face.hi.i << " || "
                  << face.lo.j << "\t" << face.hi.j << " || "
                  << face.lo.k << "\t" << face.hi.k << "\n"
                  << std::flush;
        // 如果都不是，先简单报错，便于调试
        throw std::runtime_error("detect_direction: cannot determine direction from node_box in build_inner_vertex_pattern");
    };

    auto int_to_direction = [&](int direction)
    {
        switch (direction)
        {
        case +1:
            return Direction::XPlus;
        case -1:
            return Direction::XMinus;
        case +2:
            return Direction::YPlus;
        case -2:
            return Direction::YMinus;
        case +3:
            return Direction::ZPlus;
        case -3:
            return Direction::ZMinus;
        default:
            throw std::runtime_error("int_to_direction: invalid direction");
        }
    };

    if (loc != StaggerLocation::Cell)
        return;

    PatternKey key{loc, nghost};
    if (inner_vertex_patterns_.count(key))
        return;

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;

    for (TOPO::VertexPatch &vp : topo_->inner_vertex_patches)
    {
        //=====================================================================
        // 获取邻居block中edge的dir1 dir2
        const int nb_b = vp.nb_block;
        // 取出 block
        const Block &blk_nb = fld_->grd->grids(nb_b);
        // 取出Edge的范围，以node为坐标
        Box3 tar_face_range = vp.nb_box_node;
        // tar的 cell 范围
        Int3 tar_blk = block_node_size(blk_nb);
        // 通过 this_box_node 判断接口在本块的方向 XMinus XPlus...
        int tar_dir1, tar_dir2, tar_dir3;
        detect_direction(tar_face_range, tar_blk, tar_dir1, tar_dir2, tar_dir3);
        // 保证dir1为neighbor block与this block接触面的法向
        if (vp.trans.perm[abs(vp.dir1) - 1] + 1 != abs(tar_dir1))
        {
            if (vp.trans.perm[abs(vp.dir1) - 1] + 1 == abs(tar_dir2))
            {
                int temp = tar_dir1;
                tar_dir1 = tar_dir2;
                tar_dir2 = temp;
            }
            else if (vp.trans.perm[abs(vp.dir1) - 1] + 1 == abs(tar_dir3))
            {
                int temp = tar_dir1;
                tar_dir1 = tar_dir3;
                tar_dir3 = temp;
            }
            else
            {
                std::cout << "Error for inner Edge Processing, the corresponding relation is broken! !\n";
                exit(-1);
            }
        }

        //=====================================================================
        // 构造HaloRegion
        HaloRegion r;
        r.this_block = vp.this_block;   // 顶点隶属块，接收 ghost
        r.neighbor_block = vp.nb_block; // 提供 inner
        r.this_rank = vp.this_rank;
        r.neighbor_rank = vp.nb_rank;

        r.trans = vp.trans; // this -> nb

        Direction this_d1 = int_to_direction(vp.dir1);
        Direction this_d2 = int_to_direction(vp.dir2);
        Direction this_d3 = int_to_direction(vp.dir3);

        Direction tar_d1 = int_to_direction(tar_dir1);
        Direction tar_d2 = int_to_direction(tar_dir2);
        Direction tar_d3 = int_to_direction(tar_dir3);

        // recv_box：this_block 顶点 ghost 角点
        r.recv_box = make_cell_vertex_ghost_box(vp.this_box_node, this_d1, this_d2, this_d3, nghost);

        // send_box：neighbor_block 对应 inner 角点, tar_d1为接触面法向，取inner，其他ghost
        r.send_box = make_cell_vertex_innerghost_box(vp.nb_box_node, tar_d1, tar_d2, tar_d3, nghost);

        pat.regions.push_back(r);
    }

    inner_vertex_patterns_[key] = std::move(pat);
}

void Halo::build_inner_edge_pattern(StaggerLocation loc, int nghost)
{
    auto block_node_size = [](const Block &blk) -> Int3
    {
        return {blk.mx, blk.my, blk.mz};
    };

    auto detect_direction = [](Box3 &face, Int3 &blk_mxyz, int &dir1, int &dir2)
    {
        dir1 = -999;
        dir2 = -999;
        // 下面的判断逻辑基于假设：
        // 接口在某一方向上厚度为 1 层节点，而且刚好贴在 0 或 Ni_node 上。

        // XMinus: i=0 这一层节点
        if (face.lo.i == 0 && face.hi.i == 1)
        {
            int direction = -1;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
            {
                dir2 = direction;
                return;
            }
        }

        // XPlus: i=Ni_node 这一层节点
        if (face.lo.i == blk_mxyz.i && face.hi.i == blk_mxyz.i + 1)
        {
            int direction = 1;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
            {
                dir2 = direction;
                return;
            }
        }

        // YMinus: j=0
        if (face.lo.j == 0 && face.hi.j == 1)
        {
            int direction = -2;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
            {
                dir2 = direction;
                return;
            }
        }

        // YPlus: j=Nj_node
        if (face.lo.j == blk_mxyz.j && face.hi.j == blk_mxyz.j + 1)
        {
            int direction = 2;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
            {
                dir2 = direction;
                return;
            }
        }

        // ZMinus: k=0
        if (face.lo.k == 0 && face.hi.k == 1)
        {
            int direction = -3;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
            {
                dir2 = direction;
                return;
            }
        }

        // ZPlus: k=Nk_node
        if (face.lo.k == blk_mxyz.k && face.hi.k == blk_mxyz.k + 1)
        {
            int direction = 3;
            if (dir1 == -999)
                dir1 = direction;
            else if (dir2 == -999)
            {
                dir2 = direction;
                return;
            }
        }

        std::cout << "block:" << blk_mxyz.i << " " << blk_mxyz.j << " " << blk_mxyz.k << "\n"
                  << face.lo.i << "\t" << face.hi.i << " || "
                  << face.lo.j << "\t" << face.hi.j << " || "
                  << face.lo.k << "\t" << face.hi.k << "\n"
                  << std::flush;
        // 如果都不是，先简单报错，便于调试
        throw std::runtime_error("detect_direction: cannot determine direction from node_box in build_inner_edge_pattern");
    };

    auto int_to_direction = [&](int direction)
    {
        switch (direction)
        {
        case +1:
            return Direction::XPlus;
        case -1:
            return Direction::XMinus;
        case +2:
            return Direction::YPlus;
        case -2:
            return Direction::YMinus;
        case +3:
            return Direction::ZPlus;
        case -3:
            return Direction::ZMinus;
        default:
            throw std::runtime_error("int_to_direction: invalid direction");
        }
    };

    // 暂时只给 Cell 场做角区，如果别的 Stagger 以后需要，再扩展
    if (loc != StaggerLocation::Cell)
        return;

    PatternKey key{loc, nghost};
    if (inner_edge_patterns_.find(key) != inner_edge_patterns_.end())
        return; // 已经建过

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;
    pat.regions.clear();

    for (TOPO::EdgePatch &ep : topo_->inner_edge_patches)
    {
        //=====================================================================
        // 获取邻居block中edge的dir1 dir2
        const int nb_b = ep.nb_block;
        // 取出 block
        const Block &blk_nb = fld_->grd->grids(nb_b);
        // 取出Edge的范围，以node为坐标
        Box3 tar_face_range = ep.nb_box_node;
        // tar的 cell 范围
        Int3 tar_blk = block_node_size(blk_nb);
        // 通过 this_box_node 判断接口在本块的方向 XMinus XPlus...
        int tar_dir1, tar_dir2;
        detect_direction(tar_face_range, tar_blk, tar_dir1, tar_dir2);
        // 保证dir1为neighbor block与this block接触面的法向
        if (ep.trans.perm[abs(ep.dir1) - 1] + 1 != abs(tar_dir1))
        {
            if (ep.trans.perm[abs(ep.dir1) - 1] + 1 == abs(tar_dir2))
            {
                // 安全检测
                int temp = tar_dir1;
                tar_dir1 = tar_dir2;
                tar_dir2 = temp;
            }
            else
            {
                std::cout << "Error for inner Edge Processing, the corresponding relation is broken! !\n";
                exit(-1);
            }
        }

        //=====================================================================
        // 构造HaloRegion
        HaloRegion r;

        // edge 隶属块 = 接收方（this_block）
        r.this_block = ep.this_block;   // 这里填 ghost
        r.neighbor_block = ep.nb_block; // 这里提供 inner
        r.this_rank = ep.this_rank;
        r.neighbor_rank = ep.nb_rank;

        // trans: 按 Halo_Type 的语义，是 this -> nb 的映射
        r.trans = ep.trans;

        Direction this_d1 = int_to_direction(ep.dir1);
        Direction this_d2 = int_to_direction(ep.dir2);

        Direction tar_d1 = int_to_direction(tar_dir1);
        Direction tar_d2 = int_to_direction(tar_dir2);

        // recv_box 在 this_block（edge 所属块）上
        // 注意，recv的两个方向都是虚网格故可直接调用
        r.recv_box = make_cell_edge_ghost_box(ep.this_box_node, this_d1, this_d2, nghost);

        // send_box 在 neighbor_block 上，从邻居的 inner strip 取数据
        // 注意tar_d1为两块交界面的法向，处理角区时应该取inner
        r.send_box = make_cell_edge_innerghost_box(ep.nb_box_node, tar_d1, tar_d2, nghost);

        pat.regions.push_back(r);
    }

    inner_edge_patterns_[key] = std::move(pat);
}
