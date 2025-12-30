#include "3_field/3_MPCNS_Halo.h"

Box3 Halo::make_cell_edge_ghost_box(Box3 &edge_node, Direction &dir1, Direction &dir2, int nghost)
{
    // edge_node 是“沿棱的一条 node strip”：[lo, hi)
    // 先把它映射到 Cell 内部 DOF 区间，然后在两个方向上扩展 ghost
    Box3 box;
    box.lo = edge_node.lo;
    // 和 make_cell_ghost_box 一样，把 node hi 转成 cell hi_int
    box.hi = {edge_node.hi.i - 1,
              edge_node.hi.j - 1,
              edge_node.hi.k - 1};
    // 注意：hi 是半开区间右端点（cell 内部是 [0,Ni)）

    // 在两个出界方向上依次改 ghost 范围
    switch_range_ghost(box, dir1, nghost);
    switch_range_ghost(box, dir2, nghost);

    return box;
}

// 发送则一定是一个为inner一个为ghost，约定，dir1为inner
Box3 Halo::make_cell_edge_innerghost_box(Box3 &edge_node, Direction &dir1, Direction &dir2, int nghost)
{
    // 和 ghost 逻辑类似，只是改成 inner strip
    Box3 box;
    box.lo = edge_node.lo;
    box.hi = {edge_node.hi.i - 1,
              edge_node.hi.j - 1,
              edge_node.hi.k - 1};

    // 先在第一个方向取 inner strip，再在第二个方向取 ghost strip
    switch_range_inner(box, dir1, nghost);
    switch_range_ghost(box, dir2, nghost);

    return box;
}

Box3 Halo::make_cell_vertex_ghost_box(Box3 &vertex_node, Direction &dir1, Direction &dir2, Direction &dir3, int nghost)
{
    // vertex_node 是一个 [i,i+1)×[j,j+1)×[k,k+1) 的小 node 立方
    // 先映射到内部 cell 上的那个顶点 cell，再沿三个方向扩展 ghost
    Box3 box;
    box.lo = vertex_node.lo;
    box.hi = {vertex_node.hi.i - 1,
              vertex_node.hi.j - 1,
              vertex_node.hi.k - 1};

    // 三个方向依次做 ghost 扩展，就得到三维角区:
    // 例如 XMinus+YMinus+ZMinus -> i∈[-g,0), j∈[-g,0), k∈[-g,0)
    switch_range_ghost(box, dir1, nghost);
    switch_range_ghost(box, dir2, nghost);
    switch_range_ghost(box, dir3, nghost);

    return box;
}

// 发送则一定是一个为inner两个为ghost，约定，dir1为inner
Box3 Halo::make_cell_vertex_innerghost_box(Box3 &vertex_node, Direction &dir1, Direction &dir2, Direction &dir3, int nghost)
{
    // 同理，三维顶点处的 inner 角区
    Box3 box;
    box.lo = vertex_node.lo;
    box.hi = {vertex_node.hi.i - 1,
              vertex_node.hi.j - 1,
              vertex_node.hi.k - 1};

    // 三个方向依次裁剪为 inner strip：
    // 例如 XMinus+YMinus+ZMinus -> i∈[0,g), j∈[0,g), k∈[0,g)
    switch_range_inner(box, dir1, nghost);
    switch_range_ghost(box, dir2, nghost);
    switch_range_ghost(box, dir3, nghost);

    return box;
}
