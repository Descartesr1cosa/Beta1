#pragma once
#include <vector>
#include <string>
#include <unordered_map>

#include "1_grid/1_MPCNS_Grid.h"   // Block
#include "3_field/1_Field_Block.h" // FieldBlock
#include "3_field/Field_Type.h"    // FieldDescriptor
#include "3_field/Coupling_Type.h"

class Field
{
public:
    Field() = default;
    ~Field() = default;

    Field(Grid *grd_, Param *par_)
    {
        grd = grd_;
        par = par_;
        set_blocks(grd);
        build_geometry();
    };

    // 注册一个物理场（记录 desc，立刻分配）
    void register_field(const FieldDescriptor &desc);

    // 分配所有 fieldid 下 block 的数据
    void allocate(int32_t fieldID);

    int field_id(std::string field_name) { return name_to_id_[field_name]; }
    int num_fields() const { return static_cast<int>(field_descs_.size()); }
    int num_blocks() const { return static_cast<int>(blocks_.size()); }

    const FieldDescriptor &descriptor(int32_t fid) const { return field_descs_[fid]; }

    // 按 ID 访问所有block
    std::vector<FieldBlock> &field(int32_t fid)
    {
        return field_blocks_[fid];
    }
    // 按 ID 访问
    FieldBlock &field(int32_t fid, int iblock)
    {
        return field_blocks_[fid][iblock];
    }
    // 按名字访问所有block
    std::vector<FieldBlock> &field(std::string name)
    {
        return field(name_to_id_.at(name));
    }
    // 按名字访问
    FieldBlock &field(std::string name, int iblock)
    {
        return field(name_to_id_.at(name), iblock);
    }

    //===================================================================================
    void build_geometry();
    //===================================================================================

    // 注册一个“物理对 -> 多通道”的耦合定义（src->dst 有向）
    void register_coupling_pair(const CouplingPairDesc &desc);
    // 查询（后面分配缓冲/Halo 会用）
    bool has_coupling_pair(const std::string &src, const std::string &dst) const;
    const CouplingPairDesc &coupling_pair(const std::string &src, const std::string &dst) const;

    void build_coupling_buffers(const TOPO::Topology &topo, int dimension);

private:
    // 存储网格指针
    void set_blocks(Grid *grd);

    // 这个 rank 上的所有 Block（只存指针，不拥有）
    std::vector<Block *> blocks_;

    // 所有场的描述
    std::vector<FieldDescriptor> field_descs_;
    std::unordered_map<std::string, int32_t> name_to_id_;

    std::unordered_map<std::string, std::vector<int>> blocks_by_name_;

    // 真正的数据：field_blocks_[fid][iblock]
    std::vector<std::vector<FieldBlock>> field_blocks_;

    // physic pair唯一表
    using PairKey = std::pair<std::string, std::string>; // (src,dst)
    std::map<PairKey, CouplingPairDesc> coupling_pairs_;

    std::map<PairKey, CouplingBuffersForPair> coupling_buffers_;

public:
    Grid *grd;
    Param *par;
};