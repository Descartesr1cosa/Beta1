#include "1_Boundary.h"
#include "0_basic/Error.h"

void MercuryBoundary::AddGroup(const BoundGroup &g)
{
    if (!halo_)
        ERROR::Abort("AddGroup: call Setup first");

    groups_[g.name] = g;

    if (g.do_halo)
    {
        for (auto &fn : g.fields)
            halo_->register_halo_field(fn, g.halo_level);
    }
}

void MercuryBoundary::RegisterPhysical_(const std::string &field, const std::string &region, PhysicalHandler h)
{
    bound_.RegisterPhysical(field, region, std::move(h));
}

void MercuryBoundary::RegisterCoupling_(const std::string &src, const std::string &dst,
                                        StaggerLocation loc,
                                        const std::string &channel_tag,
                                        const std::string &dst_field,
                                        CouplingHandler h)
{
    bound_.RegisterCoupling(src, dst, loc, channel_tag, dst_field, std::move(h));
}
