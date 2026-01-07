#include "3_field/Coupling_Type.h"
#include "2_topology/2_MPCNS_Topology.h"
#include "6_boundary/Boundary.h"

// ------------------------------------------------------------
// Registry: Coupling
// ------------------------------------------------------------
void BoundaryCore::RegisterCoupling(const std::string &src,
                                    const std::string &dst,
                                    StaggerLocation loc,
                                    const std::string &channel_tag,
                                    const std::string &dst_field_name,
                                    BOUND::CouplingHandler h)
{
    BOUND::CouplingKey k;
    k.src = src;
    k.dst = dst;
    k.location = loc;
    k.channel_tag = channel_tag;
    k.dst_field_name = dst_field_name;
    cpl_reg_[k] = std::move(h);
}

void BoundaryCore::RegisterCoupling(const std::string &src,
                                    const std::string &dst,
                                    StaggerLocation loc,
                                    BOUND::CouplingHandler h)
{
    // wildcard: channel_tag="" dst_field_name=""
    RegisterCoupling(src, dst, loc, "", "", std::move(h));
}

void BoundaryCore::RegisterCoupling(const std::string &src,
                                    const std::string &dst,
                                    BOUND::CouplingHandler h)
{
    // pair-level default：loc 这里无法用枚举表示 “Any”，所以用一种常见做法：
    // 约定：location=StaggerLocation::Cell 且 channel_tag="" dst_field_name="" 作为 pair default
    // 你也可以改成：只提供 per-location default，不提供 pair default。
    RegisterCoupling(src, dst, StaggerLocation::Cell, "", "", std::move(h));
}

// ------------------------------------------------------------
// Apply Coupling
// ------------------------------------------------------------
void BoundaryCore::ApplyCouplingPair(const std::string &src, const std::string &dst)
{
    // 需要 coupling buffers 已经 build 且 halo 已经把数据搬到 buf.data 里
    auto &bs = fld_->coupling_buffers(src, dst);
    const auto &channels = bs.desc.channels;

    for (int cid = 0; cid < (int)channels.size(); ++cid)
    {
        const auto &ch = channels[cid];
        const std::string &tag = ch.tag;
        const StaggerLocation loc = ch.location;

        // 默认假设：dst field name == channel tag
        // 如果你未来要 tag!=dst_field，就在 RegisterCoupling 时用 dst_field_name 覆盖，并在这里按你的规则映射。
        const std::string dst_field = tag;

        const int fid_dst = fld_->field_id(dst_field);
        (void)fid_dst; // 防止未使用告警（如果你想加 assert 可用它）

        // Resolve handler（不绑定）
        auto h = ResolveCoupling(src, dst, loc, tag, dst_field);

        auto apply_one_list = [&](std::vector<CouplingBufferBlock> &lst)
        {
            for (auto &buf : lst)
            {
                if (!buf.allocated)
                    continue;

                const int ib = buf.this_block;
                FieldBlock &Udst = fld_->field(fid_dst, ib);

                if (h)
                    h(Udst, fld_, buf, src, dst, tag);
                else
                    DefaultCouplingCopy(Udst, fld_, buf, src, dst, tag);
            }
        };

        // face
        apply_one_list(bs.inner_face[cid]);
        apply_one_list(bs.parallel_face[cid]);

        // edge
        apply_one_list(bs.inner_edge[cid]);
        apply_one_list(bs.parallel_edge[cid]);

        // vertex
        apply_one_list(bs.inner_vertex[cid]);
        apply_one_list(bs.parallel_vertex[cid]);
    }
}

BOUND::CouplingHandler BoundaryCore::ResolveCoupling(const std::string &src,
                                                     const std::string &dst,
                                                     StaggerLocation loc,
                                                     const std::string &channel_tag,
                                                     const std::string &dst_field_name) const
{
    // 优先级（建议）：
    // 1) (src,dst,loc,channel,dst_field)
    // 2) (src,dst,loc,channel,"")
    // 3) (src,dst,loc,"","")
    //
    // 可选：pair-level default（如果你想用 RegisterCoupling(src,dst,h)）
    // 这里为了不引入 AnyLocation，我们给一个弱约定 fallback：
    // 4) (src,dst,Cell,"","") 作为 pair default
    auto find_one = [&](StaggerLocation L,
                        const std::string &ch,
                        const std::string &df) -> BOUND::CouplingHandler
    {
        BOUND::CouplingKey k;
        k.src = src;
        k.dst = dst;
        k.location = L;
        k.channel_tag = ch;
        k.dst_field_name = df;

        auto it = cpl_reg_.find(k);
        if (it != cpl_reg_.end())
            return it->second;
        return nullptr;
    };

    if (auto h = find_one(loc, channel_tag, dst_field_name))
        return h;
    if (auto h = find_one(loc, channel_tag, ""))
        return h;
    if (auto h = find_one(loc, "", ""))
        return h;
    if (auto h = find_one(StaggerLocation::Cell, "", ""))
        return h; // pair default（弱约定）
    return nullptr;
}

void BoundaryCore::DefaultCouplingCopy(FieldBlock &Udst, Field * /*fld*/,
                                       CouplingBufferBlock &buf,
                                       const std::string & /*src*/,
                                       const std::string & /*dst*/,
                                       const std::string & /*channel_tag*/)
{
    // 把 buf.box 区域的数据写入 Udst 同样的索引范围
    const Box3 &b = buf.box;
    const int ncomp = buf.ncomp;

    for (int i = b.lo.i; i < b.hi.i; ++i)
        for (int j = b.lo.j; j < b.hi.j; ++j)
            for (int k = b.lo.k; k < b.hi.k; ++k)
            {
                for (int m = 0; m < ncomp; ++m)
                {
                    Udst(i, j, k, m) = buf(i, j, k, m);
                }
            }
}