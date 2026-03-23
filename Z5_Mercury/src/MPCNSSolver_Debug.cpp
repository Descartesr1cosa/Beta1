#include "MercurySolver.h"

#include <mpi.h>
#include <set>
#include <tuple>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>

void MercurySolver::DebugPrintFieldByTopo_(int query_rank,
                                           int blk_num,
                                           int i0, int j0, int k0,
                                           std::string &field_name)
{
    struct QueryItem
    {
        int rank;
        int blk;
        int i, j, k;

        int src_rank;
        int src_blk;
        int src_i, src_j, src_k;

        int patch_tag; // 0=SRC, 1=face-inner, 2=face-parallel, 3=edge-inner, 4=edge-parallel, 5=vertex-inner, 6=vertex-parallel
    };

    auto get_comp = [](const Int3 &x, int d) -> int
    {
        if (d == 0)
            return x.i;
        if (d == 1)
            return x.j;
        return x.k;
    };

    auto set_comp = [](Int3 &x, int d, int v)
    {
        if (d == 0)
            x.i = v;
        else if (d == 1)
            x.j = v;
        else
            x.k = v;
    };

    auto in_box = [](const Box3 &box, const Int3 &p) -> bool
    {
        return (box.lo.i <= p.i && p.i < box.hi.i) &&
               (box.lo.j <= p.j && p.j < box.hi.j) &&
               (box.lo.k <= p.k && p.k < box.hi.k);
    };

    auto loc_name = [](StaggerLocation loc) -> const char *
    {
        switch (loc)
        {
        case StaggerLocation::Cell:
            return "Cell";
        case StaggerLocation::Node:
            return "Node";
        case StaggerLocation::FaceXi:
            return "FaceXi";
        case StaggerLocation::FaceEt:
            return "FaceEt";
        case StaggerLocation::FaceZe:
            return "FaceZe";
        case StaggerLocation::EdgeXi:
            return "EdgeXi";
        case StaggerLocation::EdgeEt:
            return "EdgeEt";
        case StaggerLocation::EdgeZe:
            return "EdgeZe";
        default:
            return "Unknown";
        }
    };

    auto patch_name = [](int tag) -> const char *
    {
        switch (tag)
        {
        case 0:
            return "SRC";
        case 1:
            return "FACE-INNER";
        case 2:
            return "FACE-PAR";
        case 3:
            return "EDGE-INNER";
        case 4:
            return "EDGE-PAR";
        case 5:
            return "VERT-INNER";
        case 6:
            return "VERT-PAR";
        default:
            return "UNKNOWN";
        }
    };

    auto delta_on_axis = [](StaggerLocation loc, int ax) -> int
    {
        switch (loc)
        {
        case StaggerLocation::Cell:
            return 1;
        case StaggerLocation::Node:
            return 0;
        case StaggerLocation::FaceXi:
            return (ax == 0 ? 0 : 1);
        case StaggerLocation::FaceEt:
            return (ax == 1 ? 0 : 1);
        case StaggerLocation::FaceZe:
            return (ax == 2 ? 0 : 1);
        case StaggerLocation::EdgeXi:
            return (ax == 0 ? 1 : 0);
        case StaggerLocation::EdgeEt:
            return (ax == 1 ? 1 : 0);
        case StaggerLocation::EdgeZe:
            return (ax == 2 ? 1 : 0);
        default:
            return 1;
        }
    };

    // 把 topo 的 node-box 转成该 field/stagger 在“本块逻辑索引空间”里的 dof box
    auto make_dof_box_from_node_box = [&](StaggerLocation loc, const Box3 &node_box) -> Box3
    {
        Box3 box;
        box.lo = node_box.lo;

        switch (loc)
        {
        case StaggerLocation::Cell:
            box.hi = {node_box.hi.i - 1, node_box.hi.j - 1, node_box.hi.k - 1};
            break;
        case StaggerLocation::Node:
            box.hi = {node_box.hi.i, node_box.hi.j, node_box.hi.k};
            break;
        case StaggerLocation::FaceXi:
            box.hi = {node_box.hi.i, node_box.hi.j - 1, node_box.hi.k - 1};
            break;
        case StaggerLocation::FaceEt:
            box.hi = {node_box.hi.i - 1, node_box.hi.j, node_box.hi.k - 1};
            break;
        case StaggerLocation::FaceZe:
            box.hi = {node_box.hi.i - 1, node_box.hi.j - 1, node_box.hi.k};
            break;
        case StaggerLocation::EdgeXi:
            box.hi = {node_box.hi.i - 1, node_box.hi.j, node_box.hi.k};
            break;
        case StaggerLocation::EdgeEt:
            box.hi = {node_box.hi.i, node_box.hi.j - 1, node_box.hi.k};
            break;
        case StaggerLocation::EdgeZe:
            box.hi = {node_box.hi.i, node_box.hi.j, node_box.hi.k - 1};
            break;
        }
        return box;
    };

    // 把某个 staggered DOF 的锚点索引映射到邻居块
    // 公式：
    //   nb[perm[a]] = sign[a] * local[a] + offset[a]
    // 对 sign[a] = -1 的情况，要额外减去该 DOF 在该轴上的 support 长度 delta[a]
    auto map_dof_index = [&](const Int3 &p_local,
                             StaggerLocation loc,
                             const TOPO::IndexTransform &tr) -> Int3
    {
        int lv[3] = {p_local.i, p_local.j, p_local.k};
        Int3 p_nb{0, 0, 0};

        for (int a = 0; a < 3; ++a)
        {
            const int b = tr.perm[a];
            const int s = tr.sign[a];
            const int off = get_comp(tr.offset, a);
            const int delta = delta_on_axis(loc, a);

            int v = s * lv[a] + off;
            if (s < 0)
                v -= delta;

            set_comp(p_nb, b, v);
        }
        return p_nb;
    };

    auto print_value = [&](const QueryItem &q)
    {
        if (!fld_->has_field(field_name))
        {
            std::cout << "[rank " << par_->GetInt("myid") << "] "
                      << "field \"" << field_name << "\" not registered\n";
            return;
        }

        const int fid = fld_->field_id(field_name);
        const auto &desc = fld_->descriptor(fid);

        if (q.blk < 0 || q.blk >= fld_->num_blocks())
        {
            std::cout << "[rank " << par_->GetInt("myid") << "] "
                      << "[" << patch_name(q.patch_tag) << "] "
                      << "invalid blk=" << q.blk
                      << " for field=" << field_name << "\n";
            return;
        }

        auto &F = fld_->field(fid, q.blk);
        if (!F.is_allocated())
        {
            std::cout << "[rank " << par_->GetInt("myid") << "] "
                      << "[" << patch_name(q.patch_tag) << "] "
                      << "field=" << field_name
                      << " blk=" << q.blk
                      << " not allocated\n";
            return;
        }

        const Int3 lo = F.get_lo();
        const Int3 hi = F.get_hi();

        if (q.i < lo.i || q.i >= hi.i ||
            q.j < lo.j || q.j >= hi.j ||
            q.k < lo.k || q.k >= hi.k)
        {
            std::cout << "[rank " << par_->GetInt("myid") << "] "
                      << "[" << patch_name(q.patch_tag) << "] "
                      << "field=" << field_name
                      << " loc=" << loc_name(desc.location)
                      << " blk=" << q.blk
                      << " idx=(" << q.i << "," << q.j << "," << q.k << ") "
                      << "out_of_range, valid=[("
                      << lo.i << "," << lo.j << "," << lo.k << "),("
                      << hi.i << "," << hi.j << "," << hi.k << "))\n";
            return;
        }

        std::ostringstream oss;
        oss << std::setprecision(16);
        oss << "[rank " << par_->GetInt("myid") << "] "
            << "[" << patch_name(q.patch_tag) << "] "
            << "field=" << field_name
            << " loc=" << loc_name(desc.location)
            << " blk=" << q.blk
            << " idx=(" << q.i << "," << q.j << "," << q.k << ")"
            << " <- src(rank=" << q.src_rank
            << ", blk=" << q.src_blk
            << ", idx=(" << q.src_i << "," << q.src_j << "," << q.src_k << "))"
            << " : ";

        for (int m = 0; m < desc.ncomp; ++m)
        {
            oss << F(q.i, q.j, q.k, m);
            if (m + 1 < desc.ncomp)
                oss << ", ";
        }
        oss << "\n";

        std::cout << oss.str();
    };

    const int myrank = par_->GetInt("myid");
    int nrank = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    if (!fld_->has_field(field_name))
    {
        if (myrank == query_rank)
        {
            std::cout << "[DebugPrintFieldByTopo_] field \"" << field_name
                      << "\" not registered\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    const int fid = fld_->field_id(field_name);
    const auto &desc = fld_->descriptor(fid);

    std::vector<QueryItem> queries;
    std::set<std::tuple<int, int, int, int, int, int>> uniq;

    auto push_query = [&](int rank, int blk, int i, int j, int k, int patch_tag)
    {
        auto key = std::make_tuple(rank, blk, i, j, k, patch_tag);
        if (!uniq.insert(key).second)
            return;

        QueryItem q;
        q.rank = rank;
        q.blk = blk;
        q.i = i;
        q.j = j;
        q.k = k;
        q.src_rank = query_rank;
        q.src_blk = blk_num;
        q.src_i = i0;
        q.src_j = j0;
        q.src_k = k0;
        q.patch_tag = patch_tag;
        queries.push_back(q);
    };

    if (myrank == query_rank)
    {
        push_query(query_rank, blk_num, i0, j0, k0, 0);

        const Int3 p0{i0, j0, k0};

        auto scan_face_list = [&](const auto &plist, int patch_tag)
        {
            for (const auto &p : plist)
            {
                if (p.this_rank != query_rank)
                    continue;
                if (p.this_block != blk_num)
                    continue;

                const Box3 dof_box = make_dof_box_from_node_box(desc.location, p.this_box_node);
                if (!in_box(dof_box, p0))
                    continue;

                const Int3 p_nb = map_dof_index(p0, desc.location, p.trans);
                push_query(p.nb_rank, p.nb_block, p_nb.i, p_nb.j, p_nb.k, patch_tag);
            }
        };

        auto scan_edge_list = [&](const auto &plist, int patch_tag)
        {
            for (const auto &p : plist)
            {
                if (p.this_rank != query_rank)
                    continue;
                if (p.this_block != blk_num)
                    continue;

                const Box3 dof_box = make_dof_box_from_node_box(desc.location, p.this_box_node);
                if (!in_box(dof_box, p0))
                    continue;

                const Int3 p_nb = map_dof_index(p0, desc.location, p.trans);
                push_query(p.nb_rank, p.nb_block, p_nb.i, p_nb.j, p_nb.k, patch_tag);
            }
        };

        auto scan_vertex_list = [&](const auto &plist, int patch_tag)
        {
            for (const auto &p : plist)
            {
                if (p.this_rank != query_rank)
                    continue;
                if (p.this_block != blk_num)
                    continue;

                const Box3 dof_box = make_dof_box_from_node_box(desc.location, p.this_box_node);
                if (!in_box(dof_box, p0))
                    continue;

                const Int3 p_nb = map_dof_index(p0, desc.location, p.trans);
                push_query(p.nb_rank, p.nb_block, p_nb.i, p_nb.j, p_nb.k, patch_tag);
            }
        };

        scan_face_list(topo_->inner_patches, 1);
        scan_face_list(topo_->parallel_patches, 2);
        scan_edge_list(topo_->inner_edge_patches, 3);
        scan_edge_list(topo_->parallel_edge_patches, 4);
        scan_vertex_list(topo_->inner_vertex_patches, 5);
        scan_vertex_list(topo_->parallel_vertex_patches, 6);

        // 物理边界这里只做提示，不做对侧输出
        for (const auto &p : topo_->physical_patches)
        {
            if (p.this_rank != query_rank)
                continue;
            if (p.this_block != blk_num)
                continue;

            const Box3 dof_box = make_dof_box_from_node_box(desc.location, p.this_box_node);
            if (!in_box(dof_box, p0))
                continue;

            std::cout << "[rank " << myrank << "] "
                      << "[PHYSICAL] "
                      << "field=" << field_name
                      << " loc=" << loc_name(desc.location)
                      << " src blk=" << blk_num
                      << " idx=(" << i0 << "," << j0 << "," << k0 << ") "
                      << " touches physical bc \"" << p.bc_name
                      << "\" dir=" << p.direction << "\n";
        }
    }

    int nquery = (myrank == query_rank) ? static_cast<int>(queries.size()) : 0;
    MPI_Bcast(&nquery, 1, MPI_INT, query_rank, MPI_COMM_WORLD);

    if (myrank != query_rank)
        queries.resize(nquery);

    if (nquery > 0)
    {
        MPI_Bcast(reinterpret_cast<void *>(queries.data()),
                  nquery * static_cast<int>(sizeof(QueryItem)),
                  MPI_BYTE,
                  query_rank,
                  MPI_COMM_WORLD);
    }

    // 按 rank 顺序串行打印，避免 stdout 打架
    for (int r = 0; r < nrank; ++r)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        if (myrank == r)
        {
            for (const auto &q : queries)
            {
                if (q.rank != myrank)
                    continue;
                print_value(q);
            }
            std::cout.flush();
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void MercurySolver::DebugPrintEdgeEquivClass_(int query_rank,
                                              int blk, int i, int j, int k, int dir) const
{
    const int myrank = par_->GetInt("myid");
    int nrank = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);

    auto print_edge = [](const TOPO::EdgeLocalID &e) -> std::string
    {
        std::ostringstream oss;
        oss << "("
            << "rank=" << e.rank
            << ", blk=" << e.gblock
            << ", i=" << e.i
            << ", j=" << e.j
            << ", k=" << e.k
            << ", dir=" << e.dir
            << ")";
        return oss.str();
    };

    auto print_node = [](const TOPO::NodeEqID &n) -> std::string
    {
        std::ostringstream oss;
        oss << "("
            << "rank=" << n.rank
            << ", blk=" << n.gblock
            << ", i=" << n.i
            << ", j=" << n.j
            << ", k=" << n.k
            << ")";
        return oss.str();
    };

    // 先由源 rank 构造查询 edge 和其 canonical key
    TOPO::EdgeLocalID query_edge{};
    TOPO::EdgeKey key{};
    int found_on_src = 0;

    if (myrank == query_rank)
    {
        query_edge = TOPO::EdgeLocalID{query_rank, blk, i, j, k, dir};

        auto it = topo_equiv_->edge2key.find(query_edge);
        if (it != topo_equiv_->edge2key.end())
        {
            key = it->second;
            found_on_src = 1;
        }
    }

    MPI_Bcast(&found_on_src, 1, MPI_INT, query_rank, MPI_COMM_WORLD);

    if (!found_on_src)
    {
        if (myrank == query_rank)
        {
            std::cout << "[DebugEdgeEquivClass] source edge not found in edge2key: "
                      << "(rank=" << query_rank
                      << ", blk=" << blk
                      << ", i=" << i
                      << ", j=" << j
                      << ", k=" << k
                      << ", dir=" << dir << ")\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    // 广播 EdgeKey
    MPI_Bcast(&key, sizeof(TOPO::EdgeKey), MPI_BYTE, query_rank, MPI_COMM_WORLD);

    if (myrank == query_rank)
    {
        std::cout << "\n===== Edge equivalence class =====\n";
        std::cout << "source = " << print_edge(query_edge) << "\n";
        std::cout << "canonical key:\n";
        std::cout << "  a = " << print_node(key.a) << "\n";
        std::cout << "  b = " << print_node(key.b) << "\n";

        auto it_owner = topo_equiv_->edge_owner.find(key);
        if (it_owner != topo_equiv_->edge_owner.end())
        {
            const auto &eo = it_owner->second;
            std::cout << "owner rep = " << print_edge(eo);

            auto it_gid = topo_equiv_->edge_owner_gid.find(eo);
            if (it_gid != topo_equiv_->edge_owner_gid.end())
                std::cout << "  gid=" << it_gid->second;

            std::cout << "\n";
        }
        else
        {
            std::cout << "owner rep = <not found>\n";
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // 各 rank 输出自己本地属于这个 equivalence class 的成员
    for (int r = 0; r < nrank; ++r)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        if (myrank != r)
            continue;

        auto it = topo_equiv_->edge_members.find(key);
        if (it != topo_equiv_->edge_members.end())
        {
            auto members = it->second;
            std::sort(members.begin(), members.end());

            std::cout << "[rank " << myrank << "] local members for this edge class:\n";

            for (const auto &e : members)
            {
                int sign = 999;
                auto it_sign = topo_equiv_->edge2sign.find(e);
                if (it_sign != topo_equiv_->edge2sign.end())
                    sign = static_cast<int>(it_sign->second);

                bool is_owner = false;
                auto it_owner_flag = topo_equiv_->edge_is_owner.find(e);
                if (it_owner_flag != topo_equiv_->edge_is_owner.end())
                    is_owner = it_owner_flag->second;

                int gid = -1;
                auto it_gid = topo_equiv_->edge_owner_gid.find(e);
                if (it_gid != topo_equiv_->edge_owner_gid.end())
                    gid = it_gid->second;

                std::cout << "  " << print_edge(e)
                          << "  sign_to_canonical=" << sign
                          << "  is_owner=" << is_owner
                          << "  gid=" << gid
                          << "\n";
            }
        }

        std::cout.flush();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}