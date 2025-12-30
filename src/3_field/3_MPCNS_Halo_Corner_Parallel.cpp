#include "3_field/3_MPCNS_Halo.h"
#include "0_basic/MPI_WRAPPER.h"

void Halo::mpi_exchange_edge_meta(
    const std::map<int, std::vector<EdgeMeta>> &meta_to_send,
    std::vector<EdgeMeta> &recv_metas)
{
    int nrank;
    PARALLEL::mpi_size(&nrank);
    int myrank;
    PARALLEL::mpi_rank(&myrank);

    // 1) 每个 rank 要发给别人的数量数组 send_counts[nrank] / recv_counts[nrank]
    std::vector<int> send_counts(nrank, 0), recv_counts(nrank, 0);
    for (auto &kv : meta_to_send)
    {
        int dest = kv.first;
        send_counts[dest] = (int)kv.second.size();
    }

    // 用 Alltoall 或者一圈 Send/Recv 换得 recv_counts
    PARALLEL::mpi_alltoall(send_counts.data(), 1, recv_counts.data(), 1);

    // 2) 计算偏移、总数
    std::vector<int> sdispls(nrank, 0), rdispls(nrank, 0);
    int total_send = 0, total_recv = 0;
    for (int r = 0; r < nrank; ++r)
    {
        sdispls[r] = total_send;
        rdispls[r] = total_recv;
        total_send += send_counts[r];
        total_recv += recv_counts[r];
    }

    std::vector<EdgeMeta> send_buf(total_send);
    recv_metas.resize(total_recv);

    // 3) 把 meta_to_send 展平到 send_buf（按 rank 顺序）
    for (int r = 0; r < nrank; ++r)
    {
        auto it = meta_to_send.find(r);
        if (it == meta_to_send.end())
            continue;
        const auto &v = it->second;
        std::copy(v.begin(), v.end(), send_buf.begin() + sdispls[r]);
    }

    // 4) Alltoallv 交换元数据（用 MPI_BYTE 发送 struct）
    // 把“struct 个数”转换成“字节数”，调用 MPI_Alltoallv
    const int sz = static_cast<int>(sizeof(EdgeMeta));

    std::vector<int> send_counts_bytes(nrank), recv_counts_bytes(nrank);
    std::vector<int> sdispls_bytes(nrank), rdispls_bytes(nrank);
    for (int r = 0; r < nrank; ++r)
    {
        send_counts_bytes[r] = send_counts[r] * sz;
        recv_counts_bytes[r] = recv_counts[r] * sz;
        sdispls_bytes[r] = sdispls[r] * sz;
        rdispls_bytes[r] = rdispls[r] * sz;
    }

    // 5) 真正交换 meta 的 Alltoallv
    MPI_Alltoallv(
        // send
        reinterpret_cast<const void *>(send_buf.data()),
        send_counts_bytes.data(),
        sdispls_bytes.data(),
        MPI_BYTE,
        // recv
        reinterpret_cast<void *>(recv_metas.data()),
        recv_counts_bytes.data(),
        rdispls_bytes.data(),
        MPI_BYTE,
        // comm
        MPI_COMM_WORLD);
}

void Halo::build_parallel_edge_pattern(StaggerLocation loc, int nghost)
{
    int myid;
    PARALLEL::mpi_rank(&myid);

    if (loc != StaggerLocation::Cell)
        return;

    PatternKey key{loc, nghost};
    if (parallel_edge_patterns_recv.count(key))
        return;

    HaloPattern pat_recv;
    pat_recv.location = loc;
    pat_recv.nghost = nghost;

    // 每个 neighbor_rank 一个 meta 列表
    std::map<int, std::vector<EdgeMeta>> meta_to_send;

    // 每个 neighbor_rank 自增 tag
    std::map<int, int> next_tag;

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

    auto map_dir_to_neighbor = [&](int dir_local, const TOPO::IndexTransform &tr, bool is_norm) -> int
    {
        int axis_local = std::abs(dir_local) - 1;
        int axis_nb = tr.perm[axis_local];
        int s_local = (dir_local > 0) ? tr.sign[axis_local] : -tr.sign[axis_local]; // sign为翻转，dir_local正负为大小号面
        int s_multi = (is_norm) ? -s_local : s_local;                               // 法向面需要额外翻转
        return s_multi * (axis_nb + 1);
    };

    for (TOPO::EdgePatch &ep : topo_->parallel_edge_patches)
    {
        //=========================
        // 1. 先搞清发送侧的方向
        //=========================
        // 本块上的 edge 方向（两个出界轴）
        int this_dir1 = ep.dir1;
        int this_dir2 = ep.dir2;

        // 映射到邻居块坐标的方向
        int tar_dir1_i = map_dir_to_neighbor(this_dir1, ep.trans, true);
        int tar_dir2_i = map_dir_to_neighbor(this_dir2, ep.trans, false);

        // 保证 tar_dir1 是「接触面的法向」，和 ep.dir1 变换一致
        if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 != std::abs(tar_dir1_i))
        {
            if (ep.trans.perm[std::abs(ep.dir1) - 1] + 1 == std::abs(tar_dir2_i))
                std::swap(tar_dir1_i, tar_dir2_i);
            else
            {
                std::cout << "Error for parallel Edge Processing, mapping broken!\n";
                std::exit(-1);
            }
        }

        Direction this_d1 = int_to_direction(ep.dir1);
        Direction this_d2 = int_to_direction(ep.dir2);
        Direction send_d1 = int_to_direction(tar_dir1_i);
        Direction send_d2 = int_to_direction(tar_dir2_i);

        //=========================
        // 2. 构造 recv 侧 HaloRegion
        //=========================
        HaloRegion r;
        r.this_block = ep.this_block;   // recv
        r.neighbor_block = ep.nb_block; // send
        r.this_rank = ep.this_rank;
        r.neighbor_rank = ep.nb_rank;
        r.trans = ep.trans; // recv -> send

        // 接收块上的 ghost edge 区域
        r.recv_box = make_cell_edge_ghost_box(ep.this_box_node, this_d1, this_d2, nghost);

        // recv 侧其实不会用 send_box，这里可以留空
        r.send_box = Box3{};

        // 给这一条 edge 分配一个本地 tag
        int nb_rank = ep.nb_rank;
        int tag = next_tag[nb_rank]++; // 默认 0 起，不重复就行

        r.send_flag = tag;
        r.recv_flag = tag;

        pat_recv.regions.push_back(r);

        //=========================
        // 3. 打一个 EdgeMeta，准备发给 send_rank
        //=========================
        EdgeMeta meta;
        meta.key = key;
        meta.recv_rank = ep.this_rank;
        meta.send_rank = ep.nb_rank;
        meta.recv_block = ep.this_block;
        meta.send_block = ep.nb_block;
        meta.edge_node_on_send = ep.nb_box_node;
        meta.dir1_send = tar_dir1_i;
        meta.dir2_send = tar_dir2_i;
        meta.trans_recv_to_send = ep.trans; // recv->send
        meta.tag = tag;

        meta_to_send[ep.nb_rank].push_back(meta);
    }

    parallel_edge_patterns_recv[key] = std::move(pat_recv);

    //=========================
    // 4. MPI 交换 EdgeMeta
    //=========================
    std::vector<EdgeMeta> recv_metas;
    mpi_exchange_edge_meta(meta_to_send, recv_metas);

    //=========================
    // 5. 根据收到的 meta 构建 parallel_edge_patterns_send
    //=========================
    HaloPattern pat_send;
    pat_send.location = loc;
    pat_send.nghost = nghost;

    auto inverse_transform = [&](const TOPO::IndexTransform &tr) -> TOPO::IndexTransform
    {
        TOPO::IndexTransform inv;
        // nb[perm[a]] = sign[a]*loc[a] + offset[a]
        // => loc[a] = sign[a]*(nb[perm[a]] - offset[a])
        // 可以推到 inv.perm / inv.sign / inv.offset
        int inv_perm[3];
        for (int a = 0; a < 3; ++a)
            inv_perm[tr.perm[a]] = a;

        for (int b = 0; b < 3; ++b)
        {
            int a = inv_perm[b]; // nb 的第 b 轴，对应原来的 local 第 a 轴

            inv.perm[b] = a;          // nb[b] -> local[a]
            inv.sign[b] = tr.sign[a]; // 系数同 sign[a]
            // offset[a] 是原来公式 nb[...] = sign[a]*local[a] + offset[a] 里的 offset[a]
            // 逆公式里变成 -sign[a]*offset[a]
            int off_a = (a == 0   ? tr.offset.i
                         : a == 1 ? tr.offset.j
                                  : tr.offset.k);
            int off_b = -tr.sign[a] * off_a;

            if (b == 0)
                inv.offset.i = off_b;
            if (b == 1)
                inv.offset.j = off_b;
            if (b == 2)
                inv.offset.k = off_b;
        }
        return inv;
    };

    for (const EdgeMeta &m : recv_metas)
    {

        if (m.send_rank != myid)
        {
            std::cout << "Fatal Error, received data is not for my process id: " << myid << "\t but data is for " << m.send_rank << std::endl;
            exit(-1);
        }

        HaloRegion r;
        r.this_block = m.send_block; // 在发送侧，本块是 send_block
        r.neighbor_block = m.recv_block;
        r.this_rank = m.send_rank;
        r.neighbor_rank = m.recv_rank;

        // 对发送侧来说，需要 trans: this(send) -> nb(recv)
        r.trans = inverse_transform(m.trans_recv_to_send);

        Direction send_d1 = int_to_direction(m.dir1_send);
        Direction send_d2 = int_to_direction(m.dir2_send);

        // send_box：发送侧 inner strip + ghost strip（按约定，dir1 为 inner）
        r.send_box = make_cell_edge_innerghost_box(
            const_cast<Box3 &>(m.edge_node_on_send), send_d1, send_d2, nghost);

        // 发送侧通常不需要用 recv_box，这里可以空着
        r.recv_box = Box3{};

        // tag 必须和 recv 侧一致
        r.send_flag = m.tag;
        r.recv_flag = m.tag;

        pat_send.regions.push_back(r);
    }

    parallel_edge_patterns_send[key] = std::move(pat_send);
}

void Halo::mpi_exchange_vertex_meta(
    const std::map<int, std::vector<VertexMeta>> &meta_to_send,
    std::vector<VertexMeta> &recv_metas)
{
    int nrank;
    PARALLEL::mpi_size(&nrank);
    int myrank;
    PARALLEL::mpi_rank(&myrank);

    // 1) 每个 rank 要发给别人的数量数组 send_counts[nrank] / recv_counts[nrank]
    std::vector<int> send_counts(nrank, 0), recv_counts(nrank, 0);
    for (auto &kv : meta_to_send)
    {
        int dest = kv.first;
        send_counts[dest] = (int)kv.second.size();
    }

    // 用 Alltoall 或者一圈 Send/Recv 换得 recv_counts
    PARALLEL::mpi_alltoall(send_counts.data(), 1, recv_counts.data(), 1);

    // 2) 计算偏移、总数
    std::vector<int> sdispls(nrank, 0), rdispls(nrank, 0);
    int total_send = 0, total_recv = 0;
    for (int r = 0; r < nrank; ++r)
    {
        sdispls[r] = total_send;
        rdispls[r] = total_recv;
        total_send += send_counts[r];
        total_recv += recv_counts[r];
    }

    std::vector<VertexMeta> send_buf(total_send);
    recv_metas.resize(total_recv);

    // 3) 把 meta_to_send 展平到 send_buf（按 rank 顺序）
    for (int r = 0; r < nrank; ++r)
    {
        auto it = meta_to_send.find(r);
        if (it == meta_to_send.end())
            continue;
        const auto &v = it->second;
        std::copy(v.begin(), v.end(), send_buf.begin() + sdispls[r]);
    }

    // 4) Alltoallv 交换元数据（用 MPI_BYTE 发送 struct）
    // 把“struct 个数”转换成“字节数”，调用 MPI_Alltoallv
    const int sz = static_cast<int>(sizeof(VertexMeta));

    std::vector<int> send_counts_bytes(nrank), recv_counts_bytes(nrank);
    std::vector<int> sdispls_bytes(nrank), rdispls_bytes(nrank);
    for (int r = 0; r < nrank; ++r)
    {
        send_counts_bytes[r] = send_counts[r] * sz;
        recv_counts_bytes[r] = recv_counts[r] * sz;
        sdispls_bytes[r] = sdispls[r] * sz;
        rdispls_bytes[r] = rdispls[r] * sz;
    }

    // 5) 真正交换 meta 的 Alltoallv
    MPI_Alltoallv(
        // send
        reinterpret_cast<const void *>(send_buf.data()),
        send_counts_bytes.data(),
        sdispls_bytes.data(),
        MPI_BYTE,
        // recv
        reinterpret_cast<void *>(recv_metas.data()),
        recv_counts_bytes.data(),
        rdispls_bytes.data(),
        MPI_BYTE,
        // comm
        MPI_COMM_WORLD);
}

void Halo::build_parallel_vertex_pattern(StaggerLocation loc, int nghost)
{
    int myid;
    PARALLEL::mpi_rank(&myid);

    if (loc != StaggerLocation::Cell)
        return;

    PatternKey key{loc, nghost};
    if (parallel_vertex_patterns_recv.count(key))
        return;

    HaloPattern pat_recv;
    pat_recv.location = loc;
    pat_recv.nghost = nghost;

    // 每个 neighbor_rank 一个 meta 列表
    std::map<int, std::vector<VertexMeta>> meta_to_send;

    // 每个 neighbor_rank 自增 tag
    std::map<int, int> next_tag;

    auto block_node_size = [](const Block &blk) -> Int3
    {
        return {blk.mx, blk.my, blk.mz};
    };

    // 和 inner_vertex 的 detect_direction 一样
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

    auto map_dir_to_neighbor = [&](int dir_local, const TOPO::IndexTransform &tr, bool is_norm) -> int
    {
        int axis_local = std::abs(dir_local) - 1;
        int axis_nb = tr.perm[axis_local];
        int s_local = (dir_local > 0) ? tr.sign[axis_local] : -tr.sign[axis_local]; // sign为翻转，dir_local正负为大小号面
        int s_multi = (is_norm) ? -s_local : s_local;                               // 法向面需要额外翻转
        return s_multi * (axis_nb + 1);
    };

    for (TOPO::VertexPatch &vp : topo_->parallel_vertex_patches)
    {
        //=========================
        // 1. 先搞清发送侧的方向
        //=========================
        // 本块上的 edge 方向（两个出界轴）
        int this_dir1 = vp.dir1;
        int this_dir2 = vp.dir2;
        int this_dir3 = vp.dir3;

        // 映射到邻居块坐标的方向
        int tar_dir1_i = map_dir_to_neighbor(this_dir1, vp.trans, true);
        int tar_dir2_i = map_dir_to_neighbor(this_dir2, vp.trans, false);
        int tar_dir3_i = map_dir_to_neighbor(this_dir3, vp.trans, false);

        // 保证 tar_dir1 是「接触面的法向」，和 ep.dir1 变换一致
        if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 != std::abs(tar_dir1_i))
        {
            if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(tar_dir2_i))
                std::swap(tar_dir1_i, tar_dir2_i);
            else if (vp.trans.perm[std::abs(vp.dir1) - 1] + 1 == std::abs(tar_dir3_i))
                std::swap(tar_dir1_i, tar_dir3_i);
            else
            {
                std::cout << "Error for parallel Vertex Processing, mapping broken!\n";
                std::exit(-1);
            }
        }

        Direction this_d1 = int_to_direction(vp.dir1);
        Direction this_d2 = int_to_direction(vp.dir2);
        Direction this_d3 = int_to_direction(vp.dir3);

        Direction send_d1 = int_to_direction(tar_dir1_i);
        Direction send_d2 = int_to_direction(tar_dir2_i);
        Direction send_d3 = int_to_direction(tar_dir3_i);

        //=========================
        // 2. 构造 recv 侧 HaloRegion
        //=========================
        HaloRegion r;
        r.this_block = vp.this_block;   // recv
        r.neighbor_block = vp.nb_block; // send
        r.this_rank = vp.this_rank;
        r.neighbor_rank = vp.nb_rank;
        r.trans = vp.trans; // recv -> send

        // 接收块上的 ghost edge 区域
        r.recv_box = make_cell_vertex_ghost_box(vp.this_box_node, this_d1, this_d2, this_d3, nghost);

        // recv 侧其实不会用 send_box，这里可以留空
        r.send_box = Box3{};

        // 给这一条 edge 分配一个本地 tag
        int nb_rank = vp.nb_rank;
        int tag = next_tag[nb_rank]++; // 默认 0 起，不重复就行

        r.send_flag = tag;
        r.recv_flag = tag;

        pat_recv.regions.push_back(r);

        //=========================
        // 3. 打一个 EdgeMeta，准备发给 send_rank
        //=========================
        VertexMeta meta;
        meta.key = key;
        meta.recv_rank = vp.this_rank;
        meta.send_rank = vp.nb_rank;
        meta.recv_block = vp.this_block;
        meta.send_block = vp.nb_block;
        meta.vertex_node_on_send = vp.nb_box_node;
        meta.dir1_send = tar_dir1_i;
        meta.dir2_send = tar_dir2_i;
        meta.dir3_send = tar_dir3_i;
        meta.trans_recv_to_send = vp.trans; // recv->send
        meta.tag = tag;

        meta_to_send[vp.nb_rank].push_back(meta);
    }

    parallel_vertex_patterns_recv[key] = std::move(pat_recv);

    //=========================
    // 4. MPI 交换 EdgeMeta
    //=========================
    std::vector<VertexMeta> recv_metas;
    mpi_exchange_vertex_meta(meta_to_send, recv_metas);

    //=========================
    // 5. 根据收到的 meta 构建 parallel_edge_patterns_send
    //=========================
    HaloPattern pat_send;
    pat_send.location = loc;
    pat_send.nghost = nghost;

    auto inverse_transform = [&](const TOPO::IndexTransform &tr) -> TOPO::IndexTransform
    {
        TOPO::IndexTransform inv;
        // nb[perm[a]] = sign[a]*loc[a] + offset[a]
        // => loc[a] = sign[a]*(nb[perm[a]] - offset[a])
        // 可以推到 inv.perm / inv.sign / inv.offset
        int inv_perm[3];
        for (int a = 0; a < 3; ++a)
            inv_perm[tr.perm[a]] = a;

        for (int b = 0; b < 3; ++b)
        {
            int a = inv_perm[b]; // nb 的第 b 轴，对应原来的 local 第 a 轴

            inv.perm[b] = a;          // nb[b] -> local[a]
            inv.sign[b] = tr.sign[a]; // 系数同 sign[a]
            // offset[a] 是原来公式 nb[...] = sign[a]*local[a] + offset[a] 里的 offset[a]
            // 逆公式里变成 -sign[a]*offset[a]
            int off_a = (a == 0   ? tr.offset.i
                         : a == 1 ? tr.offset.j
                                  : tr.offset.k);
            int off_b = -tr.sign[a] * off_a;

            if (b == 0)
                inv.offset.i = off_b;
            if (b == 1)
                inv.offset.j = off_b;
            if (b == 2)
                inv.offset.k = off_b;
        }
        return inv;
    };

    for (const VertexMeta &m : recv_metas)
    {

        if (m.send_rank != myid)
        {
            std::cout << "Fatal Error, received data is not for my process id: " << myid << "\t but data is for " << m.send_rank << std::endl;
            exit(-1);
        }

        HaloRegion r;
        r.this_block = m.send_block; // 在发送侧，本块是 send_block
        r.neighbor_block = m.recv_block;
        r.this_rank = m.send_rank;
        r.neighbor_rank = m.recv_rank;

        // 对发送侧来说，需要 trans: this(send) -> nb(recv)
        r.trans = inverse_transform(m.trans_recv_to_send);

        Direction send_d1 = int_to_direction(m.dir1_send);
        Direction send_d2 = int_to_direction(m.dir2_send);
        Direction send_d3 = int_to_direction(m.dir3_send);

        // send_box：发送侧 inner strip + ghost strip（按约定，dir1 为 inner）
        r.send_box = make_cell_vertex_innerghost_box(
            const_cast<Box3 &>(m.vertex_node_on_send), send_d1, send_d2, send_d3, nghost);

        // 发送侧通常不需要用 recv_box，这里可以空着
        r.recv_box = Box3{};

        // tag 必须和 recv 侧一致
        r.send_flag = m.tag;
        r.recv_flag = m.tag;

        pat_send.regions.push_back(r);
    }

    parallel_vertex_patterns_send[key] = std::move(pat_send);
}