#include "ibtk/DebuggingUtilities.h"

#include "LS/LSFindCellVolume.h"
#include "LS/utility_functions.h"

#include "ibamr/app_namespaces.h"

namespace LS
{
LSFindCellVolume::LSFindCellVolume(std::string object_name, Pointer<PatchHierarchy<NDIM>> hierarchy)
    : d_object_name(std::move(object_name)), d_hierarchy(hierarchy)
{
    // intentionally blank
    return;
} // Constructor
} // namespace LS
