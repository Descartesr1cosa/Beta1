#include "4_halo/1_MPCNS_Halo.h"
#include <set>

void Halo::register_halo_field(const std::string &field_name, HaloLevel level)
{
    // 1) 防止拼错名字导致 silent 插入
    if (!fld_->has_field(field_name))
    {
        std::cout << "[Halo] ERROR: register_halo_field: field_name = "
                  << field_name << " not registered in Field.\n";
        std::exit(-1);
    }

    // 2) 升级策略：同名多次注册取更高等级
    auto it = halo_registry_.find(field_name);
    if (it == halo_registry_.end())
    {
        halo_registry_.emplace(field_name, level);
    }
    else
    {
        int old_lv = static_cast<int>(it->second);
        int new_lv = static_cast<int>(level);
        it->second = static_cast<HaloLevel>(std::max(old_lv, new_lv));
    }
}

void Halo::build_registered_patterns()
{
    // 清理旧 pattern（可重复 build）
    inner_patterns_.clear();
    parallel_patterns_.clear();
    inner_edge_patterns_.clear();
    parallel_edge_patterns_send.clear();
    parallel_edge_patterns_recv.clear();
    inner_vertex_patterns_.clear();
    parallel_vertex_patterns_send.clear();
    parallel_vertex_patterns_recv.clear();

    // coupling parallel corner patterns
    coupling_parallel_edge_patterns_send.clear();
    coupling_parallel_edge_patterns_recv.clear();
    coupling_parallel_vertex_patterns_send.clear();
    coupling_parallel_vertex_patterns_recv.clear();

    // 维度
    const int dim = fld_->grd->dimension;

    // 1) Build patterns per field. Material-restricted fields can share
    //    location/nghost with full-domain fields, but their active block set
    //    is different.
    for (const auto &kv : halo_registry_)
    {
        const std::string &fname = kv.first;
        HaloLevel lv = kv.second;

        int fid = fld_->field_id(fname); // 这里要求 field_id 不产生副作用；最好配合 has_field
        const auto &desc = fld_->descriptor(fid);

        build_inner_1DCorner_pattern(fname, desc.location, desc.nghost);
        build_parallel_1DCorner_pattern(fname, desc.location, desc.nghost);

        if (dim >= 2 && static_cast<int>(lv) >= static_cast<int>(HaloLevel::Edge))
        {
            build_inner_2DCorner_pattern(fname, desc.location, desc.nghost);
            build_parallel_2DCorner_pattern(fname, desc.location, desc.nghost);
        }

        if (dim >= 3 && static_cast<int>(lv) >= static_cast<int>(HaloLevel::Vertex))
        {
            build_inner_3DCorner_pattern(fname, desc.location, desc.nghost);
            build_parallel_3DCorner_pattern(fname, desc.location, desc.nghost);
        }
    }

    // 2) Coupling parallel corner patterns (directed src -> dst)
    //    Build once here to avoid lazy rebuild during frequent coupling exchanges.
    const auto &cpairs = fld_->coupling_pairs();
    if (!cpairs.empty())
    {
        // small helpers: check whether any matching coupling patches exist on this rank
        auto has_parallel_coupling_edge = [&](const std::string &src, const std::string &dst) -> bool
        {
            for (const auto &ep : topo_->parallel_edge_patches)
                if (ep.is_coupling && ep.nb_block_name == src && ep.this_block_name == dst)
                    return true;
            return false;
        };
        auto has_parallel_coupling_vertex = [&](const std::string &src, const std::string &dst) -> bool
        {
            for (const auto &vp : topo_->parallel_vertex_patches)
                if (vp.is_coupling && vp.nb_block_name == src && vp.this_block_name == dst)
                    return true;
            return false;
        };

        for (const auto &kv : cpairs)
        {
            const CouplingPairDesc &pd = kv.second;
            const std::string &src = pd.pair.src;
            const std::string &dst = pd.pair.dst;

            std::set<std::pair<StaggerLocation, int>> ckeys;
            for (const auto &ch : pd.channels)
                ckeys.insert({ch.location, ch.nghost});

            if (dim >= 2 && has_parallel_coupling_edge(src, dst))
            {
                for (const auto &k : ckeys)
                    build_coupling_parallel_2DCorner_pattern(src, dst, k.first, k.second);
            }

            if (dim >= 3 && has_parallel_coupling_vertex(src, dst))
            {
                for (const auto &k : ckeys)
                    build_coupling_parallel_3DCorner_pattern(src, dst, k.first, k.second);
            }
        }
    }
}
