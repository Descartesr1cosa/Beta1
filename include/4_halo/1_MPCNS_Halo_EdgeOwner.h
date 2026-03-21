#pragma once
#include <cstdint>
#include <vector>

#include "2_topology/2_MPCNS_Topology_Equiv.h"
#include "3_field/2_MPCNS_Field.h"
#include "4_halo/Halo_EdgeOwner_Type.h"

namespace HALO_OWNER
{
    // ============================================================
    // edge-owner sync pattern
    // ============================================================

    struct EdgeOwnerLocalAliasItem
    {
        TOPO::EdgeLocalID owner; // local owner rep
        TOPO::EdgeLocalID rep;   // local non-owner rep
        int8_t sign;             // for 1-form; ignored by vec copy
    };

    struct EdgeOwnerSendItem
    {
        int tar_id;              // target rank id
        TOPO::EdgeLocalID owner; // local owner rep
    };

    struct EdgeOwnerRecvItem
    {
        int tar_id;            // remote owner rank id
        TOPO::EdgeLocalID rep; // local non-owner rep
        int8_t sign;           // for 1-form; ignored by vec copy
    };

    struct EdgeOwnerSyncPattern
    {
        std::vector<EdgeOwnerLocalAliasItem> local_alias;

        std::vector<EdgeOwnerSendItem> send_items;
        std::vector<EdgeOwnerRecvItem> recv_items;

        std::vector<int> send_counts;
        std::vector<int> recv_counts;
        std::vector<int> send_displs;
        std::vector<int> recv_displs;

        void clear()
        {
            local_alias.clear();
            send_items.clear();
            recv_items.clear();

            send_counts.clear();
            recv_counts.clear();
            send_displs.clear();
            recv_displs.clear();
        }
    };

    // ============================================================
    // build
    // ============================================================

    void build_edge_owner_sync_pattern(
        const TOPO::TopologyEquiv &equiv,
        EdgeOwnerSyncPattern &pattern);

    // ============================================================
    // runtime sync
    // ============================================================

    // Treat every component on edge as a 1-form component:
    // rep(comp) = sign * owner(comp)
    void sync_edge_1form(
        Field &fld,
        IdTriplet &field_id,
        const EdgeOwnerSyncPattern &pattern);

    // Treat every component on edge as a Cartesian vector component:
    // rep(comp) = owner(comp)
    void sync_edge_vec(
        Field &fld,
        IdTriplet &field_id,
        const EdgeOwnerSyncPattern &pattern);

    void pack_owner_edge_1form_local(
        Field &fld,
        const IdTriplet &field_id,
        const TOPO::TopologyEquiv &equiv,
        std::vector<double> &buf_local);

    void unpack_owner_edge_1form_local(
        const std::vector<double> &buf_local,
        Field &fld,
        const IdTriplet &field_id,
        const TOPO::TopologyEquiv &equiv,
        const EdgeOwnerSyncPattern &pattern);

} // namespace HALO_OWNER
