#include "3_field/3_MPCNS_Halo.h"

void Halo::build_pattern()
{
    for (int fid = 0; fid < fld_->num_fields(); fid++)
    {
        auto loc = fld_->descriptor(fid).location;
        auto nghost = fld_->descriptor(fid).nghost;

        // 1. face 级
        build_inner_pattern(loc, nghost);
        build_parallel_pattern(loc, nghost);
    }

    int U_fid = fld_->field_id("U_");
    auto loc = fld_->descriptor(U_fid).location;
    auto nghost = fld_->descriptor(U_fid).nghost;
    // 2. edge 级
    build_inner_edge_pattern(loc, nghost);
    build_parallel_edge_pattern(loc, nghost);

    // 3. vertex 级
    build_inner_vertex_pattern(loc, nghost);
    build_parallel_vertex_pattern(loc, nghost);
}

void Halo::build_inner_pattern(StaggerLocation loc, int nghost)
{

    auto block_node_size = [](const Block &blk) -> Int3
    {
        return {blk.mx, blk.my, blk.mz};
    };

    auto detect_direction = [](Box3 &face, Int3 &blk_mxyz) -> Direction
    {
        // 下面的判断逻辑基于假设：
        // 接口在某一方向上厚度为 1 层节点，而且刚好贴在 0 或 Ni_node 上。

        // XMinus: i=0 这一层节点
        if (face.lo.i == 0 && face.hi.i == 1)
            return Direction::XMinus;

        // XPlus: i=Ni_node 这一层节点
        if (face.lo.i == blk_mxyz.i && face.hi.i == blk_mxyz.i + 1)
            return Direction::XPlus;

        // YMinus: j=0
        if (face.lo.j == 0 && face.hi.j == 1)
            return Direction::YMinus;

        // YPlus: j=Nj_node
        if (face.lo.j == blk_mxyz.j && face.hi.j == blk_mxyz.j + 1)
            return Direction::YPlus;

        // ZMinus: k=0
        if (face.lo.k == 0 && face.hi.k == 1)
            return Direction::ZMinus;

        // ZPlus: k=Nk_node
        if (face.lo.k == blk_mxyz.k && face.hi.k == blk_mxyz.k + 1)
            return Direction::ZPlus;

        // 如果都不是，先简单报错，便于调试
        throw std::runtime_error("detect_direction: cannot determine direction from node_box in build_inner_pattern");
    };

    PatternKey key{loc, nghost};
    // 已经 build 过了，直接返回，避免重复构造
    if (inner_patterns_.find(key) != inner_patterns_.end())
        return;

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;
    pat.regions.clear();

    // 遍历所有 inner_patches（同一 rank 的块-块接口）
    for (const TOPO::InterfacePatch &patch : topo_->inner_patches)
    {
        const int this_b = patch.this_block;
        const int nb_b = patch.nb_block;

        // 取出两个 block
        const Block &blk_this = fld_->grd->grids(this_b); // blocks_ 是 Field 里的
        const Block &blk_nb = fld_->grd->grids(nb_b);

        // 取出面的范围，以node为坐标
        Box3 my_face_range = patch.this_box_node;
        Box3 tar_face_range = patch.nb_box_node;

        // 各自的 cell 范围
        Int3 my_blk = block_node_size(blk_this);
        Int3 tar_blk = block_node_size(blk_nb);

        // 通过 this_box_node 判断接口在本块的方向 XMinus XPlus...
        Direction dir_my = detect_direction(my_face_range, my_blk);
        Direction dir_tar = detect_direction(tar_face_range, tar_blk);

        // --------- region 1: this_block 的 inner 到 nb_block 的 ghost  ---------
        HaloRegion r1;
        r1.this_block = this_b;   // 发送方
        r1.neighbor_block = nb_b; // 接收方
        r1.this_rank = patch.this_rank;
        r1.neighbor_rank = patch.nb_rank;
        r1.trans = patch.trans; // 从 this_block 逻辑坐标 -> neighbor_block 逻辑坐标

        // ---------处理范围 ---------
        switch (loc)
        {
        case StaggerLocation::Cell:
            // my_block 的inner区域 --> tar_block 的ghost 区域 ：
            r1.send_box = make_cell_inner_box(my_face_range, dir_my, nghost);
            r1.recv_box = make_cell_ghost_box(tar_face_range, dir_tar, nghost);
            break;
        case StaggerLocation::FaceXi:
            // my_block 的inner区域 --> tar_block 的ghost 区域 ：
            r1.send_box = make_facex_inner_box(my_face_range, dir_my, nghost);
            r1.recv_box = make_facex_ghost_box(tar_face_range, dir_tar, nghost);
            break;
        case StaggerLocation::FaceEt:
            // my_block 的inner区域 --> tar_block 的ghost 区域 ：
            r1.send_box = make_facey_inner_box(my_face_range, dir_my, nghost);
            r1.recv_box = make_facey_ghost_box(tar_face_range, dir_tar, nghost);
            break;
        case StaggerLocation::FaceZe:
            // my_block 的inner区域 --> tar_block 的ghost 区域 ：
            r1.send_box = make_facez_inner_box(my_face_range, dir_my, nghost);
            r1.recv_box = make_facez_ghost_box(tar_face_range, dir_tar, nghost);
            break;
        case StaggerLocation::EdgeXi:
            // my_block 的inner区域 --> tar_block 的ghost 区域 ：
            r1.send_box = make_edgex_inner_box(my_face_range, dir_my, nghost);
            r1.recv_box = make_edgex_ghost_box(tar_face_range, dir_tar, nghost);
            break;
        case StaggerLocation::EdgeEt:
            // my_block 的inner区域 --> tar_block 的ghost 区域 ：
            r1.send_box = make_edgey_inner_box(my_face_range, dir_my, nghost);
            r1.recv_box = make_edgey_ghost_box(tar_face_range, dir_tar, nghost);
            break;
        case StaggerLocation::EdgeZe:
            // my_block 的inner区域 --> tar_block 的ghost 区域 ：
            r1.send_box = make_edgez_inner_box(my_face_range, dir_my, nghost);
            r1.recv_box = make_edgez_ghost_box(tar_face_range, dir_tar, nghost);
            break;
        default:
            break;
        }

        // ---------加入regions ---------
        pat.regions.push_back(r1);
    }
    // 存储供exchange_inner使用
    inner_patterns_[key] = std::move(pat);
}

void Halo::build_parallel_pattern(StaggerLocation loc, int nghost)
{
    // 和 build_inner_pattern 里一致的小工具
    auto block_node_size = [](const Block &blk) -> Int3
    {
        return {blk.mx, blk.my, blk.mz};
    };

    auto detect_direction = [](Box3 &face, Int3 &blk_mxyz) -> Direction
    {
        // 接口在某一方向上厚度为 1 层节点，而且刚好贴在 0 或 Ni_node 上。

        // XMinus: i=0
        if (face.lo.i == 0 && face.hi.i == 1)
            return Direction::XMinus;

        // XPlus: i=Ni_node
        if (face.lo.i == blk_mxyz.i && face.hi.i == blk_mxyz.i + 1)
            return Direction::XPlus;

        // YMinus: j=0
        if (face.lo.j == 0 && face.hi.j == 1)
            return Direction::YMinus;

        // YPlus: j=Nj_node
        if (face.lo.j == blk_mxyz.j && face.hi.j == blk_mxyz.j + 1)
            return Direction::YPlus;

        // ZMinus: k=0
        if (face.lo.k == 0 && face.hi.k == 1)
            return Direction::ZMinus;

        // ZPlus: k=Nk_node
        if (face.lo.k == blk_mxyz.k && face.hi.k == blk_mxyz.k + 1)
            return Direction::ZPlus;

        throw std::runtime_error("detect_direction: cannot determine direction from node_box in build_parallel_pattern");
    };

    PatternKey key{loc, nghost};
    if (parallel_patterns_.find(key) != parallel_patterns_.end())
        return; // 已经建过，直接返回

    HaloPattern pat;
    pat.location = loc;
    pat.nghost = nghost;
    pat.regions.clear();

    // 遍历所有 parallel_patches（跨 rank 的接口）
    for (const TOPO::InterfacePatch &patch : topo_->parallel_patches)
    {
        const int this_b = patch.this_block;

        // 只知道本 rank 上的 block，邻居 block 不在本 rank 上
        const Block &blk_this = fld_->grd->grids(this_b);

        // 本块 interface 的 node 区域
        Box3 my_face_range = patch.this_box_node;
        Int3 my_blk = block_node_size(blk_this);

        // 判定在本块上的方向
        Direction dir_my = detect_direction(my_face_range, my_blk);

        HaloRegion r;
        r.this_block = this_b;             // 本块（既是发送方也是接收方，从 MPI 角度看）
        r.neighbor_block = patch.nb_block; // 对面的 block index，Parallel 情况下一般=-1，仅做记录
        r.this_rank = patch.this_rank;
        r.neighbor_rank = patch.nb_rank; // 真实的邻居 rank
        r.trans = patch.trans;           // 从 this_block 逻辑坐标 -> neighbor_block 逻辑坐标
        r.send_flag = patch.send_flag;
        r.recv_flag = patch.recv_flag;

        // 对 parallel 来说，send_box / recv_box 都在本块的 FieldBlock 空间里：
        //   recv_box: 本块 ghost 区域（MPI 接收后，往这里 unpack）
        //   send_box: 本块 inner strip（打包后，MPI 发送给 neighbor_rank）
        switch (loc)
        {
        case StaggerLocation::Cell:
            r.recv_box = make_cell_ghost_box(my_face_range, dir_my, nghost);
            r.send_box = make_cell_inner_box(my_face_range, dir_my, nghost);
            break;

        case StaggerLocation::FaceXi:
            r.recv_box = make_facex_ghost_box(my_face_range, dir_my, nghost);
            r.send_box = make_facex_inner_box(my_face_range, dir_my, nghost);
            break;

        case StaggerLocation::FaceEt:
            r.recv_box = make_facey_ghost_box(my_face_range, dir_my, nghost);
            r.send_box = make_facey_inner_box(my_face_range, dir_my, nghost);
            break;

        case StaggerLocation::FaceZe:
            r.recv_box = make_facez_ghost_box(my_face_range, dir_my, nghost);
            r.send_box = make_facez_inner_box(my_face_range, dir_my, nghost);
            break;

        case StaggerLocation::EdgeXi:
            r.recv_box = make_edgex_ghost_box(my_face_range, dir_my, nghost);
            r.send_box = make_edgex_inner_box(my_face_range, dir_my, nghost);
            break;

        case StaggerLocation::EdgeEt:
            r.recv_box = make_edgey_ghost_box(my_face_range, dir_my, nghost);
            r.send_box = make_edgey_inner_box(my_face_range, dir_my, nghost);
            break;

        case StaggerLocation::EdgeZe:
            r.recv_box = make_edgez_ghost_box(my_face_range, dir_my, nghost);
            r.send_box = make_edgez_inner_box(my_face_range, dir_my, nghost);
            break;

        default:
            // 目前暂时不支持 Node / Edge 之类的话可以先跳过
            continue;
        }

        pat.regions.push_back(r);
    }

    // 存起来，供 exchange_parallel 使用
    parallel_patterns_[key] = std::move(pat);
}

void Halo::switch_range_ghost(Box3 &box, Direction &dir, int nghost)
{
    // 这里输入的Box3 &face为特定StaggerLocation的坐标
    // Cell [0,mx-1] = [0,mx)
    switch (dir)
    {
    case Direction::XMinus:
        // 左 ghost：i ∈ [-g .. -1] = [-g,0)
        box.lo.i = -nghost;
        box.hi.i = 0;
        break;

    case Direction::XPlus:
        // 右 ghost：i ∈ [Ni .. Ni+g) = [hi_int, hi_int+g)
        box.lo.i = box.hi.i;
        box.hi.i = box.hi.i + nghost;
        break;

    case Direction::YMinus:
        box.lo.j = -nghost;
        box.hi.j = 0;
        break;

    case Direction::YPlus:
        box.lo.j = box.hi.j;
        box.hi.j = box.hi.j + nghost;
        break;

    case Direction::ZMinus:
        box.lo.k = -nghost;
        box.hi.k = 0;
        break;

    case Direction::ZPlus:
        box.lo.k = box.hi.k;
        box.hi.k = box.hi.k + nghost;
        break;
    }
}

void Halo::switch_range_inner(Box3 &box, Direction &dir, int nghost)
{
    switch (dir)
    {
    case Direction::XMinus:
        // 左 inner：i ∈ [0 .. g-1] = [0,g)
        box.lo.i = 0;
        box.hi.i = nghost;
        break;

    case Direction::XPlus:
        // 右 inner strip：i ∈ [Ni-g .. Ni) = [hi_int-g, hi_int)
        box.lo.i = box.hi.i - nghost;
        box.hi.i = box.hi.i;
        break;

    case Direction::YMinus:
        box.lo.j = 0;
        box.hi.j = nghost;
        break;

    case Direction::YPlus:
        box.lo.j = box.hi.j - nghost;
        box.hi.j = box.hi.j;
        break;

    case Direction::ZMinus:
        box.lo.k = 0;
        box.hi.k = nghost;
        break;

    case Direction::ZPlus:
        box.lo.k = box.hi.k - nghost;
        box.hi.k = box.hi.k;
        break;
    }
}

Box3 Halo::make_cell_ghost_box(Box3 &face, Direction &dir, int nghost)
{
    // face是基于node的范围[lo,hi)
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j - 1, face.hi.k - 1};
    // Cell 内部 i=[0..Ni-1], j=[0..Nj-1], k=[0..Nk-1]

    switch_range_ghost(box, dir, nghost);
    return box;
}

Box3 Halo::make_cell_inner_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j - 1, face.hi.k - 1};
    // Cell 内部 i=[0..Ni-1], j=[0..Nj-1], k=[0..Nk-1]
    switch_range_inner(box, dir, nghost);
    return box;
}

Box3 Halo::make_facex_ghost_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i, face.hi.j - 1, face.hi.k - 1};
    switch_range_ghost(box, dir, nghost);
    return box;
}

Box3 Halo::make_facex_inner_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i, face.hi.j - 1, face.hi.k - 1};
    switch_range_inner(box, dir, nghost);
    return box;
}

Box3 Halo::make_facey_ghost_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j, face.hi.k - 1};
    switch_range_ghost(box, dir, nghost);
    return box;
}

Box3 Halo::make_facey_inner_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j, face.hi.k - 1};
    switch_range_inner(box, dir, nghost);
    return box;
}

Box3 Halo::make_facez_ghost_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j - 1, face.hi.k};
    switch_range_ghost(box, dir, nghost);
    return box;
}

Box3 Halo::make_facez_inner_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j - 1, face.hi.k};
    switch_range_inner(box, dir, nghost);
    return box;
}

Box3 Halo::make_edgex_ghost_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j, face.hi.k};
    switch_range_ghost(box, dir, nghost);
    return box;
}

Box3 Halo::make_edgex_inner_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i - 1, face.hi.j, face.hi.k};
    switch_range_inner(box, dir, nghost);
    return box;
}

Box3 Halo::make_edgey_ghost_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i, face.hi.j - 1, face.hi.k};
    switch_range_ghost(box, dir, nghost);
    return box;
}

Box3 Halo::make_edgey_inner_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i, face.hi.j - 1, face.hi.k};
    switch_range_inner(box, dir, nghost);
    return box;
}

Box3 Halo::make_edgez_ghost_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i, face.hi.j, face.hi.k - 1};
    switch_range_ghost(box, dir, nghost);
    return box;
}

Box3 Halo::make_edgez_inner_box(Box3 &face, Direction &dir, int nghost)
{
    // 先设成“内部 DOF 范围”，再按方向改某一维
    Box3 box;
    box.lo = face.lo;
    box.hi = {face.hi.i, face.hi.j, face.hi.k - 1};
    switch_range_inner(box, dir, nghost);
    return box;
}
