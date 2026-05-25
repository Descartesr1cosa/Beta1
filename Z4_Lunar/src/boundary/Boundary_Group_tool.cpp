#include "1_Boundary.h"
#include "0_basic/Error.h"

void LunarBoundary::AddGroup(const BoundGroup &g)
{
    if (!halo_)
        ERROR::Abort("AddGroup: call Setup first");

    groups_[g.name] = g;
}

void LunarBoundary::RegisterPhysical_(const std::string &field, const std::string &region, PhysicalHandler h)
{
    bound_.RegisterPhysical(field, region, std::move(h));
}
