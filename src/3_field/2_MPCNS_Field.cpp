// 2_Fields.cpp
#include "3_field/2_MPCNS_Field.h"

void Field::set_blocks(Grid *grd)
{
    blocks_.resize(grd->nblock);
    for (int i = 0; i < grd->nblock; i++)
        blocks_[i] = &(grd->grids(i));

    // 如果之前已经注册了场，就重新分配一次, for safety
    if (!field_descs_.empty())
    {
        allocate_all();
    }
}

void Field::register_field(const FieldDescriptor &desc)
{
    int32_t fid = static_cast<int32_t>(field_descs_.size());
    field_descs_.push_back(desc);
    name_to_id_[desc.name] = fid;

    // 如果已经有 blocks_，可以只为这个 field 分配一遍；这里简单起见统一走 allocate_all。
    if (!blocks_.empty())
    {
        allocate(fid);
    }

    return;
}

void Field::allocate_all()
{
    const int nf = static_cast<int>(field_descs_.size());
    const int nb = static_cast<int>(blocks_.size());

    field_blocks_.clear();
    field_blocks_.resize(nf);
    for (int fid = 0; fid < nf; ++fid)
    {
        field_blocks_[fid].resize(nb);
        for (int b = 0; b < nb; ++b)
        {
            field_blocks_[fid][b].allocate(*blocks_[b], field_descs_[fid]);
        }
    }
}

void Field::allocate(int32_t fieldID)
{
    const int nb = static_cast<int>(blocks_.size());

    std::vector<FieldBlock> tmp;
    tmp.resize(nb);
    for (int b = 0; b < nb; ++b)
        tmp[b].allocate(*blocks_[b], field_descs_[fieldID]);

    field_blocks_.push_back(tmp);
}