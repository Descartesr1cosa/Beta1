#pragma once
#include <vector>
#include <string>
#include <unordered_map>

#include "1_grid/1_MPCNS_Grid.h"   // Block
#include "3_field/1_Field_Block.h" // FieldBlock
#include "3_field/Field_Type.h"    // FieldDescriptor

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

    // 分配所有 field × block 的数据
    void allocate_all();

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

private:
    // 存储网格指针
    void set_blocks(Grid *grd);

    // 这个 rank 上的所有 Block（只存指针，不拥有）
    std::vector<Block *> blocks_;

    // 所有场的描述
    std::vector<FieldDescriptor> field_descs_;
    std::unordered_map<std::string, int32_t> name_to_id_;

    // 真正的数据：field_blocks_[fid][iblock]
    std::vector<std::vector<FieldBlock>> field_blocks_;

public:
    Grid *grd;
    Param *par;
    // // 用一批 Block 指针初始化（一个 rank 上的所有本地块）
    // explicit Field(const std::vector<Block *> &blocks)
    // {
    //     reset_blocks(blocks);
    // }
    // // 如果你以后想换一批 Block（比如重分块），可以重置
    // void reset_blocks(const std::vector<Block *> &blocks);
    // // 遍历某个 block 上所有场：f(fid, desc, fb)
    // template <typename Func>
    // void for_each_field_on_block(int iblock, Func f)
    // {
    //     const int nf = num_fields();
    //     for (int32_t fid = 0; fid < nf; ++fid)
    //     {
    //         f(fid, field_descs_[fid], field_blocks_[fid][iblock]);
    //     }
    // }
};