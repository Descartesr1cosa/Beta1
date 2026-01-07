#include "6_boundary/Boundary.h"

// ------------------------------------------------------------
// Geometry helpers (same style as coupling buffer build)
// ------------------------------------------------------------
Int3 BoundaryCore::LocDelta(StaggerLocation loc)
{
    switch (loc)
    {
    case StaggerLocation::Cell:
        return {1, 1, 1};
    case StaggerLocation::Node:
        return {0, 0, 0};
    case StaggerLocation::FaceXi:
        return {0, 1, 1};
    case StaggerLocation::FaceEt:
        return {1, 0, 1};
    case StaggerLocation::FaceZe:
        return {1, 1, 0};
    case StaggerLocation::EdgeXi:
        return {1, 0, 0};
    case StaggerLocation::EdgeEt:
        return {0, 1, 0};
    case StaggerLocation::EdgeZe:
        return {0, 0, 1};
    default:
        return {0, 0, 0};
    }
}

Int3 BoundaryCore::LocInnerHi(const Block &blk, StaggerLocation loc)
{
    // blk.mx/my/mz 是 cell counts (Ni,Nj,Nk)
    const Int3 nodes = {blk.mx + 1, blk.my + 1, blk.mz + 1};
    const Int3 d = LocDelta(loc);
    return {nodes.i - d.i, nodes.j - d.j, nodes.k - d.k}; // half-open [0,hi)
}

void BoundaryCore::ConvertTangent(int lo_n, int hi_n, int delta, int &lo, int &hi)
{
    lo = lo_n;
    hi = (delta == 0) ? hi_n : (hi_n - 1);
}

Box3 BoundaryCore::MakeFaceGhostSlabBox(const Block &blk,
                                        StaggerLocation loc,
                                        const Box3 &face_node_box,
                                        int dir_code,
                                        int nghost)
{
    const int ax = std::abs(dir_code); // 1/2/3
    const int sgn = (dir_code > 0) ? +1 : -1;

    const Int3 hi_in = LocInnerHi(blk, loc);
    const Int3 d = LocDelta(loc);

    int t1, t2;
    if (ax == 1)
    {
        t1 = 2;
        t2 = 3;
    }
    else if (ax == 2)
    {
        t1 = 1;
        t2 = 3;
    }
    else
    {
        t1 = 1;
        t2 = 2;
    }

    Box3 b{};

    // normal ghost range
    if (ax == 1)
    {
        b.lo.i = (sgn < 0) ? -nghost : hi_in.i;
        b.hi.i = (sgn < 0) ? 0 : hi_in.i + nghost;
    }
    if (ax == 2)
    {
        b.lo.j = (sgn < 0) ? -nghost : hi_in.j;
        b.hi.j = (sgn < 0) ? 0 : hi_in.j + nghost;
    }
    if (ax == 3)
    {
        b.lo.k = (sgn < 0) ? -nghost : hi_in.k;
        b.hi.k = (sgn < 0) ? 0 : hi_in.k + nghost;
    }

    // tangential from node box -> loc box
    auto set_tangent = [&](int t)
    {
        int lo, hi;
        if (t == 1)
            ConvertTangent(face_node_box.lo.i, face_node_box.hi.i, d.i, lo, hi);
        if (t == 2)
            ConvertTangent(face_node_box.lo.j, face_node_box.hi.j, d.j, lo, hi);
        if (t == 3)
            ConvertTangent(face_node_box.lo.k, face_node_box.hi.k, d.k, lo, hi);

        if (t == 1)
        {
            b.lo.i = lo;
            b.hi.i = hi;
        }
        if (t == 2)
        {
            b.lo.j = lo;
            b.hi.j = hi;
        }
        if (t == 3)
        {
            b.lo.k = lo;
            b.hi.k = hi;
        }
    };
    set_tangent(t1);
    set_tangent(t2);

    return b;
}
